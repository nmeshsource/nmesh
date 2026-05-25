#!/usr/bin/env python3
"""Convert Nmesh SDMF text metadata to XDMF metadata.

This utility reads an SDMF ``.txt`` metadata file and writes an
XDMF ``.xmf`` file that points at the same binary data files.
It does not read, write, or convert the binary ``.bin`` files.
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path


B_HEAD_XMF = """<?xml version="1.0" encoding="utf-8"?>
<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.1">
  <Domain>
"""

E_HEAD_XMF = """  </Domain>
</Xdmf>
"""

B_TEMPORAL_XMF = """    <Grid CollectionType="Temporal" GridType="Collection" Name="TCollection">
      <Geometry Type="None"/>
      <Topology Dimensions="0" Type="NoTopology"/>
"""

E_TEMPORAL_XMF = """    </Grid>
"""

B_SPATIAL_XMF = """      <Grid CollectionType="Spatial" GridType="Collection" Name="SCollection">
        <Time Value="{time:.9f}"/>
        <Geometry Type="None"/>
        <Topology Dimensions="0" Type="NoTopology"/>
"""

E_SPATIAL_XMF = """      </Grid>
"""

B_E_GRID_XMF = """        <Grid Name="{name}">
          <Time Value="{time:.9f}"/>
          <Geometry Type="XYZ">
            <DataItem DataType="Float" Dimensions="{point_count} 3" Format="Binary" Seek="{xyz_seek}" Precision="4">
              {xyz_filename}
            </DataItem>
          </Geometry>
          <Topology Dimensions="{n2} {n1} {n0}" Type="3DSMesh"/>
          <Attribute Center="Node" Name="{attribute_name}" Type="Scalar">
            <DataItem DataType="Float" Dimensions="{n2} {n1} {n0}" Format="Binary" Seek="{value_seek}" Precision="4">
              {value_filename}
            </DataItem>
          </Attribute>
        </Grid>
"""


@dataclass(frozen=True)
class Element:
    name: str
    n0: int
    n1: int
    n2: int
    xyz_seek: int
    value_seek: int


@dataclass(frozen=True)
class TimeStep:
    time: float
    elements: tuple[Element, ...]


def parse_sdmf_names(path: Path) -> tuple[str, str]:
    if path.suffix != ".txt":
        raise ValueError(f"input file extension must be .txt, got {path.name}")

    without_txt = path.with_suffix("")
    suffix = without_txt.suffix.lstrip(".")
    attribute_name = without_txt.with_suffix("").name
    if not suffix or not attribute_name:
        raise ValueError(f"could not infer attribute name and suffix from {path.name}")

    return attribute_name, suffix


def parse_time(line: str) -> float | None:
    match = re.match(r'#\s+"time\s*=\s*([^"]+)"', line)
    if match is None:
        return None
    return float(match.group(1))


def parse_sdmf(path: Path) -> tuple[str, str, tuple[TimeStep, ...]]:
    attribute_name, suffix = parse_sdmf_names(path)
    timesteps: list[TimeStep] = []
    current_time: float | None = None
    current_elements: list[Element] = []

    with path.open("r", encoding="utf-8") as stream:
        header = stream.readline()
        if not header.startswith("# sdmf:"):
            raise ValueError(f"sdmf header missing in {path}")
        if "binarydata: float" not in header:
            raise ValueError("only SDMF binarydata: float is supported")
        if "TopologyType: 3DSMesh" not in header:
            raise ValueError("only SDMF TopologyType: 3DSMesh is supported")
        if "AttributeCenter: Node" not in header:
            raise ValueError("only SDMF AttributeCenter: Node is supported")

        for line_number, line in enumerate(stream, start=2):
            stripped = line.strip()
            if not stripped:
                continue

            if stripped.startswith("#"):
                time = parse_time(stripped)
                if time is None:
                    continue
                if current_time is not None:
                    timesteps.append(TimeStep(current_time, tuple(current_elements)))
                current_time = time
                current_elements = []
                continue

            if current_time is None:
                raise ValueError(f"element record before first time at line {line_number}")

            parts = stripped.split()
            if len(parts) != 6:
                raise ValueError(f"bad SDMF element record at line {line_number}: {stripped}")
            try:
                element = Element(
                    name=parts[0],
                    n0=int(parts[1]),
                    n1=int(parts[2]),
                    n2=int(parts[3]),
                    xyz_seek=int(parts[4]),
                    value_seek=int(parts[5]),
                )
            except ValueError as exc:
                raise ValueError(f"bad integer in SDMF element record at line {line_number}") from exc
            current_elements.append(element)

    if current_time is not None:
        timesteps.append(TimeStep(current_time, tuple(current_elements)))

    if not timesteps:
        raise ValueError(f"no SDMF timesteps found in {path}")

    return attribute_name, suffix, tuple(timesteps)


def sdmf_to_xmf_text(path: Path) -> str:
    attribute_name, suffix, timesteps = parse_sdmf(path)
    value_filename = f"{attribute_name}.{suffix}.bin"
    xyz_filename = f"xyz.{suffix}.bin"

    parts: list[str] = [B_HEAD_XMF, B_TEMPORAL_XMF]
    for timestep in timesteps:
        parts.append(B_SPATIAL_XMF.format(time=timestep.time))
        for element in timestep.elements:
            parts.append(
                B_E_GRID_XMF.format(
                    name=element.name,
                    time=timestep.time,
                    point_count=element.n0 * element.n1 * element.n2,
                    xyz_seek=element.xyz_seek,
                    xyz_filename=xyz_filename,
                    n0=element.n0,
                    n1=element.n1,
                    n2=element.n2,
                    attribute_name=attribute_name,
                    value_seek=element.value_seek,
                    value_filename=value_filename,
                )
            )
        parts.append(E_SPATIAL_XMF)
    parts.extend((E_TEMPORAL_XMF, E_HEAD_XMF))
    return "".join(parts)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert Nmesh SDMF .txt metadata to XDMF .xmf metadata."
    )
    parser.add_argument("txt", type=Path, help="input SDMF .txt file")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="output .xmf path; defaults to replacing .txt with .xmf",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    input_path = args.txt.resolve()
    output_path = args.output.resolve() if args.output is not None else input_path.with_suffix(".xmf")

    try:
        text = sdmf_to_xmf_text(input_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(text, encoding="utf-8")
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    print(f"wrote {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Convert Nmesh XDMF metadata to SDMF text metadata.

This converts only the metadata file. It does not read, write, or transform
the binary .bin files referenced by the XDMF/SDMF metadata.
"""

from __future__ import annotations

import argparse
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


HEADER = (
    "# sdmf:   binarydata: float   TopologyType: 3DSMesh   AttributeCenter: Node\n"
    "# elm n[0] n[1] n[2] xyzseek varseek\n"
)


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
    time_label: str
    elements: tuple[Element, ...]


def strip_namespace(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def children_named(element: ET.Element, name: str) -> Iterable[ET.Element]:
    for child in list(element):
        if strip_namespace(child.tag) == name:
            yield child


def first_child_named(element: ET.Element, name: str) -> ET.Element | None:
    return next(children_named(element, name), None)


def parse_dimensions(value: str) -> tuple[int, ...]:
    return tuple(int(part) for part in value.split())


def parse_time(element: ET.Element, fallback: float) -> tuple[float, str]:
    time = first_child_named(element, "Time")
    if time is None:
        label = format(fallback, ".9f")
        return fallback, label
    label = time.attrib["Value"]
    return float(label), label


def sidecar_time_labels(path: Path) -> list[str] | None:
    sidecar = path.with_suffix(".txt")
    if not sidecar.exists():
        return None

    labels: list[str] = []
    with sidecar.open("r", encoding="utf-8") as stream:
        for line in stream:
            time = parse_sdmf_time(line.strip())
            if time is not None:
                labels.append(time)
    return labels


def parse_data_item(parent: ET.Element, grid_name: str, item_name: str) -> ET.Element:
    item = first_child_named(parent, "DataItem")
    if item is None:
        raise ValueError(f"grid {grid_name} has no {item_name} DataItem")
    if item.attrib.get("DataType", "Float") != "Float":
        raise ValueError(f"grid {grid_name} has unsupported {item_name} DataType")
    if item.attrib.get("Format") != "Binary":
        raise ValueError(f"grid {grid_name} has unsupported {item_name} Format")
    if item.attrib.get("Precision", "4") != "4":
        raise ValueError(f"grid {grid_name} has unsupported {item_name} Precision")
    return item


def parse_xdmf(path: Path) -> tuple[tuple[TimeStep, ...], str | None]:
    tree = ET.parse(path)
    root = tree.getroot()
    domain = first_child_named(root, "Domain")
    if domain is None:
        raise ValueError("XDMF file has no Domain")
    temporal_grid = first_child_named(domain, "Grid")
    if temporal_grid is None:
        raise ValueError("XDMF file has no temporal Grid")

    timesteps: list[TimeStep] = []
    attribute_name: str | None = None

    for spatial_grid in children_named(temporal_grid, "Grid"):
        if spatial_grid.attrib.get("GridType") != "Collection":
            continue

        time, time_label = parse_time(spatial_grid, float(len(timesteps)))
        elements: list[Element] = []
        for grid in children_named(spatial_grid, "Grid"):
            grid_name = grid.attrib.get("Name", str(len(elements)))
            topology = first_child_named(grid, "Topology")
            geometry = first_child_named(grid, "Geometry")
            if topology is None or geometry is None:
                continue
            if topology.attrib.get("Type") != "3DSMesh":
                raise ValueError(f"grid {grid_name} has unsupported Topology Type")
            if geometry.attrib.get("Type") != "XYZ":
                raise ValueError(f"grid {grid_name} has unsupported Geometry Type")

            geometry_item = parse_data_item(geometry, grid_name, "geometry")

            attributes = list(children_named(grid, "Attribute"))
            node_attributes = [
                attribute
                for attribute in attributes
                if attribute.attrib.get("Center") == "Node"
                and attribute.attrib.get("Type", "Scalar") == "Scalar"
            ]
            if len(node_attributes) != 1:
                raise ValueError(f"grid {grid_name} must have exactly one scalar node Attribute")
            attribute = node_attributes[0]
            name = attribute.attrib.get("Name")
            if name is None:
                raise ValueError(f"grid {grid_name} has an Attribute without Name")
            if attribute_name is None:
                attribute_name = name
            elif attribute_name != name:
                raise ValueError("all Attributes must have the same Name")

            attribute_item = parse_data_item(attribute, grid_name, "attribute")
            n2, n1, n0 = parse_dimensions(topology.attrib["Dimensions"])
            elements.append(
                Element(
                    name=grid_name,
                    n0=n0,
                    n1=n1,
                    n2=n2,
                    xyz_seek=int(geometry_item.attrib.get("Seek", "0")),
                    value_seek=int(attribute_item.attrib.get("Seek", "0")),
                )
            )

        if elements:
            timesteps.append(TimeStep(time=time, time_label=time_label, elements=tuple(elements)))

    if not timesteps:
        raise ValueError("XDMF file has no 3DSMesh timesteps")
    return tuple(timesteps), attribute_name


def format_time(time: float) -> str:
    return format(time, ".16g")


def parse_sdmf_time(line: str) -> str | None:
    if not line.startswith('# "time = ') or not line.endswith('"'):
        return None
    return line[len('# "time = '):-1]


def xdmf_to_sdmf_text(path: Path) -> str:
    timesteps, _ = parse_xdmf(path)
    labels = sidecar_time_labels(path)
    if labels is not None:
        if len(labels) != len(timesteps):
            raise ValueError(
                f"sidecar {path.with_suffix('.txt')} has {len(labels)} times, "
                f"but XDMF has {len(timesteps)} timesteps"
            )
        time_labels = labels
    else:
        time_labels = [format_time(timestep.time) for timestep in timesteps]

    parts = [HEADER]
    for timestep, time_label in zip(timesteps, time_labels):
        parts.append(f'\n# "time = {time_label}"\n')
        for element in timestep.elements:
            parts.append(
                f"{element.name}\t{element.n0} {element.n1} {element.n2}\t"
                f"{element.xyz_seek} {element.value_seek}\n"
            )
    return "".join(parts)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert Nmesh XDMF .xmf metadata to SDMF .txt metadata."
    )
    parser.add_argument("xmf", type=Path, help="input XDMF .xmf file")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="output .txt path; defaults to replacing .xmf with .txt",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    input_path = args.xmf.resolve()
    output_path = args.output.resolve() if args.output is not None else input_path.with_suffix(".txt")

    try:
        text = xdmf_to_sdmf_text(input_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(text, encoding="utf-8")
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    print(f"wrote {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

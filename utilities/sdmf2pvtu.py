#!/usr/bin/env python3
"""Convert SDMF 3DSMesh binary output to VTK PVTU/VTU files.

This utility reads the SDMF text files produced by Nmesh, where each
time slice is a collection of 3DSMesh blocks with explicit XYZ coordinates and
node-centered scalar data in raw float32 binary files.

Copyright (C) 2026 Wolfgang Tichy
"""

from __future__ import annotations

import argparse
import os
import re
import struct
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO
from xml.sax.saxutils import escape


VTK_VERTEX = 1
VTK_LINE = 3
VTK_QUAD = 9
VTK_HEXAHEDRON = 12
DEFAULT_VTU_FILES_PER_PVTU = 8


@dataclass(frozen=True)
class DataItem:
    filename: str
    dimensions: tuple[int, ...]
    precision: int
    seek: int


@dataclass(frozen=True)
class Piece:
    name: str
    topology_dimensions: tuple[int, int, int]
    geometry: DataItem
    attribute_name: str
    attribute_components: int
    attribute: DataItem


@dataclass(frozen=True)
class TimeStep:
    index: int
    value: str
    pieces: tuple[Piece, ...]


def parse_dimensions(value: str) -> tuple[int, ...]:
    return tuple(int(part) for part in value.split())


def point_count(dimensions: tuple[int, int, int]) -> int:
    nx, ny, nz = dimensions
    return nx * ny * nz


def data_item_count(dimensions: tuple[int, ...]) -> int:
    count = 1
    for dimension in dimensions:
        count *= dimension
    return count


def attribute_components(
    attribute_name: str,
    attribute_type: str,
    data_item: DataItem,
    dimensions: tuple[int, int, int],
) -> int:
    item_count = data_item_count(data_item.dimensions)
    nodes = point_count(dimensions)
    if item_count % nodes != 0:
        raise ValueError(
            f"attribute {attribute_name} dimensions {data_item.dimensions} "
            f"do not match topology {dimensions}"
        )

    components = item_count // nodes
    declared_type = attribute_type
    expected_components = {
        "Scalar": 1,
        "Vector": 3,
        "Tensor": 9,
        "Tensor6": 6,
        "Matrix": components,
    }.get(declared_type)
    if expected_components is not None and expected_components != components:
        raise ValueError(
            f"attribute {attribute_name} Type={declared_type!r} "
            f"implies {expected_components} components but dimensions {data_item.dimensions} "
            f"contain {components}"
        )
    return components


def parse_header(line: str, path: Path) -> tuple[int, str, str]:
    if not line.startswith("# sdmf:"):
        raise ValueError(f"sdmf header missing in {path}")

    binary_match = re.search(r"\bbinarydata:\s+(\S+)", line)
    topology_match = re.search(r"\bTopologyType:\s+(\S+)", line)
    center_match = re.search(r"\bAttributeCenter:\s+(\S+)", line)
    if binary_match is None or topology_match is None or center_match is None:
        raise ValueError(f"incomplete sdmf header in {path}")

    binary_type = binary_match.group(1)
    if binary_type == "float":
        precision = 4
    elif binary_type == "double":
        precision = 8
    else:
        raise ValueError(f"unsupported SDMF binarydata type {binary_type!r}")

    return precision, topology_match.group(1), center_match.group(1)


def parse_time(line: str) -> str | None:
    match = re.match(r'#\s+"time\s*=\s*([^"]+)"', line)
    if match is None:
        return None
    return f"{float(match.group(1)):.9f}"


def sdmf_names(path: Path) -> tuple[str, str, str]:
    if path.suffix != ".txt":
        raise ValueError(f"input file extension must be .txt, got {path.name}")

    without_txt = path.with_suffix("")
    suffix = without_txt.suffix.lstrip(".")
    if not suffix:
        raise ValueError(f"input file {path.name} has no SDMF output suffix")

    attribute_name = without_txt.with_suffix("").name
    if not attribute_name:
        raise ValueError(f"could not determine attribute name from {path.name}")

    return attribute_name, suffix, without_txt.name


def parse_sdmf(path: Path, attribute_name: str | None) -> tuple[str, tuple[TimeStep, ...]]:
    discovered_attribute_name, suffix, _ = sdmf_names(path)
    if attribute_name is not None and attribute_name != discovered_attribute_name:
        raise ValueError(f"SDMF file contains attribute {discovered_attribute_name!r}, not {attribute_name!r}")

    value_filename = f"{discovered_attribute_name}.{suffix}.bin"
    geometry_filename = f"xyz.{suffix}.bin"
    timesteps: list[TimeStep] = []
    current_time: str | None = None
    current_pieces: list[Piece] = []

    with path.open("r", encoding="utf-8") as stream:
        header = stream.readline()
        precision, topology_type, attribute_center = parse_header(header, path)
        if topology_type != "3DSMesh":
            raise ValueError(f"unsupported topology {topology_type!r}")
        if attribute_center != "Node":
            raise ValueError(f"unsupported attribute center {attribute_center!r}")

        for line_number, line in enumerate(stream, start=2):
            stripped = line.strip()
            if not stripped:
                continue

            if stripped.startswith("#"):
                time_value = parse_time(stripped)
                if time_value is None:
                    continue
                if current_time is not None:
                    if not current_pieces:
                        raise ValueError(f"time {current_time} has no SDMF element records")
                    timesteps.append(
                        TimeStep(len(timesteps), current_time, tuple(current_pieces))
                    )
                current_time = time_value
                current_pieces = []
                continue

            if current_time is None:
                raise ValueError(f"element record before first time at line {line_number}")

            parts = stripped.split()
            if len(parts) != 6:
                raise ValueError(f"bad SDMF element record at line {line_number}: {stripped}")

            name = parts[0].rstrip("_")
            try:
                n0, n1, n2 = (int(parts[1]), int(parts[2]), int(parts[3]))
                xyz_seek, value_seek = (int(parts[4]), int(parts[5]))
            except ValueError as exc:
                raise ValueError(f"bad integer in SDMF element record at line {line_number}") from exc

            topology_dimensions = (n0, n1, n2)
            nodes = point_count(topology_dimensions)
            geometry = DataItem(
                filename=geometry_filename,
                dimensions=(nodes, 3),
                precision=precision,
                seek=xyz_seek,
            )
            attribute = DataItem(
                filename=value_filename,
                dimensions=topology_dimensions,
                precision=precision,
                seek=value_seek,
            )
            current_pieces.append(
                Piece(
                    name=name,
                    topology_dimensions=topology_dimensions,
                    geometry=geometry,
                    attribute_name=discovered_attribute_name,
                    attribute_components=attribute_components(
                        discovered_attribute_name,
                        "Scalar",
                        attribute,
                        topology_dimensions,
                    ),
                    attribute=attribute,
                )
            )

    if current_time is not None:
        if not current_pieces:
            raise ValueError(f"time {current_time} has no SDMF element records")
        timesteps.append(TimeStep(len(timesteps), current_time, tuple(current_pieces)))

    if not timesteps:
        raise ValueError("no SDMF timesteps found")

    return discovered_attribute_name, tuple(timesteps)


def parse_sdmf_file(path: Path, attribute_name: str | None) -> tuple[str, tuple[TimeStep, ...]]:
    return parse_sdmf(path, attribute_name)


class BinaryItemReader:
    def __init__(self, base_dir: Path):
        self.base_dir = base_dir
        self.streams: dict[str, BinaryIO] = {}

    def __enter__(self) -> "BinaryItemReader":
        return self

    def __exit__(self, exc_type: object, exc_value: object, traceback: object) -> None:
        self.close()

    def close(self) -> None:
        for stream in self.streams.values():
            stream.close()
        self.streams.clear()

    def read(self, item: DataItem) -> bytes:
        path = self.base_dir / item.filename
        size = data_item_count(item.dimensions) * item.precision
        stream = self.streams.get(item.filename)
        if stream is None:
            stream = path.open("rb")
            self.streams[item.filename] = stream

        stream.seek(item.seek)
        data = stream.read(size)

        if len(data) != size:
            raise ValueError(f"{path} ended while reading {size} bytes at offset {item.seek}")
        return data


def make_cells(dimensions: tuple[int, int, int], point_offset: int = 0) -> tuple[bytes, bytes, bytes, int]:
    nx, ny, nz = dimensions
    if nx < 1 or ny < 1 or nz < 1:
        raise ValueError(f"3DSMesh dimensions must all be at least 1, got {dimensions}")

    connectivity = bytearray()
    offsets = bytearray()
    types = bytearray()
    offset = 0

    def point_id(i: int, j: int, k: int) -> int:
        return point_offset + (k * ny + j) * nx + i

    active_axes = [axis for axis, dimension in enumerate(dimensions) if dimension > 1]
    axis_sizes = (nx, ny, nz)

    if len(active_axes) == 3:
        for k in range(nz - 1):
            for j in range(ny - 1):
                for i in range(nx - 1):
                    ids = (
                        point_id(i, j, k),
                        point_id(i + 1, j, k),
                        point_id(i + 1, j + 1, k),
                        point_id(i, j + 1, k),
                        point_id(i, j, k + 1),
                        point_id(i + 1, j, k + 1),
                        point_id(i + 1, j + 1, k + 1),
                        point_id(i, j + 1, k + 1),
                    )
                    connectivity.extend(struct.pack("<8i", *ids))
                    offset += 8
                    offsets.extend(struct.pack("<i", offset))
                    types.extend(struct.pack("<B", VTK_HEXAHEDRON))
    elif len(active_axes) == 2:
        a, b = active_axes
        coord = [0, 0, 0]
        for ib in range(axis_sizes[b] - 1):
            for ia in range(axis_sizes[a] - 1):
                ids = []
                for da, db in ((0, 0), (1, 0), (1, 1), (0, 1)):
                    coord[a] = ia + da
                    coord[b] = ib + db
                    ids.append(point_id(*coord))
                connectivity.extend(struct.pack("<4i", *ids))
                offset += 4
                offsets.extend(struct.pack("<i", offset))
                types.extend(struct.pack("<B", VTK_QUAD))
    elif len(active_axes) == 1:
        (a,) = active_axes
        coord = [0, 0, 0]
        for ia in range(axis_sizes[a] - 1):
            coord[a] = ia
            first = point_id(*coord)
            coord[a] = ia + 1
            second = point_id(*coord)
            connectivity.extend(struct.pack("<2i", first, second))
            offset += 2
            offsets.extend(struct.pack("<i", offset))
            types.extend(struct.pack("<B", VTK_LINE))
    else:
        connectivity.extend(struct.pack("<i", point_id(0, 0, 0)))
        offsets.extend(struct.pack("<i", 1))
        types.extend(struct.pack("<B", VTK_VERTEX))

    return bytes(connectivity), bytes(offsets), bytes(types), len(types)


def appended_array(offset: int, data: bytes) -> tuple[int, bytes]:
    return offset + 4 + len(data), struct.pack("<I", len(data)) + data


def write_vtu(
    path: Path,
    points: bytes,
    values: bytes,
    cells: tuple[bytes, bytes, bytes],
    point_count: int,
    cell_count: int,
    attribute_name: str,
    attribute_components: int,
    precision: int,
) -> None:
    vtk_float = "Float32" if precision == 4 else "Float64"
    blocks: list[bytes] = []
    offset = 0

    value_offset = offset
    offset, block = appended_array(offset, values)
    blocks.append(block)

    point_offset = offset
    offset, block = appended_array(offset, points)
    blocks.append(block)

    connectivity_offset = offset
    offset, block = appended_array(offset, cells[0])
    blocks.append(block)

    offsets_offset = offset
    offset, block = appended_array(offset, cells[1])
    blocks.append(block)

    types_offset = offset
    _, block = appended_array(offset, cells[2])
    blocks.append(block)

    xml = f"""<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" header_type="UInt32">
  <UnstructuredGrid>
    <Piece NumberOfPoints="{point_count}" NumberOfCells="{cell_count}">
      <PointData Scalars="{escape(attribute_name)}">
        <DataArray type="{vtk_float}" Name="{escape(attribute_name)}" NumberOfComponents="{attribute_components}" format="appended" offset="{value_offset}"/>
      </PointData>
      <CellData/>
      <Points>
        <DataArray type="{vtk_float}" NumberOfComponents="3" format="appended" offset="{point_offset}"/>
      </Points>
      <Cells>
        <DataArray type="Int32" Name="connectivity" format="appended" offset="{connectivity_offset}"/>
        <DataArray type="Int32" Name="offsets" format="appended" offset="{offsets_offset}"/>
        <DataArray type="UInt8" Name="types" format="appended" offset="{types_offset}"/>
      </Cells>
    </Piece>
  </UnstructuredGrid>
  <AppendedData encoding="raw">
_""".encode("utf-8")
    xml_tail = b"\n  </AppendedData>\n</VTKFile>\n"

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as stream:
        stream.write(xml)
        for block in blocks:
            stream.write(block)
        stream.write(xml_tail)


def write_pvtu(
    path: Path,
    piece_paths: list[Path],
    attribute_name: str,
    attribute_components: int,
    precision: int,
) -> None:
    vtk_float = "Float32" if precision == 4 else "Float64"
    pieces = "\n".join(
        f'      <Piece Source="{escape(os.path.relpath(piece, path.parent))}"/>'
        for piece in piece_paths
    )
    text = f"""<?xml version="1.0"?>
<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian" header_type="UInt32">
  <PUnstructuredGrid GhostLevel="0">
    <PPointData Scalars="{escape(attribute_name)}">
      <PDataArray type="{vtk_float}" Name="{escape(attribute_name)}" NumberOfComponents="{attribute_components}"/>
    </PPointData>
    <PCellData/>
    <PPoints>
      <PDataArray type="{vtk_float}" NumberOfComponents="3"/>
    </PPoints>
{pieces}
  </PUnstructuredGrid>
</VTKFile>
"""
    path.write_text(text, encoding="utf-8")


def write_pvd(path: Path, datasets: list[tuple[str, Path]]) -> None:
    rows = "\n".join(
        f'    <DataSet timestep="{escape(time)}" group="" part="0" file="{escape(os.path.relpath(pvtu, path.parent))}"/>'
        for time, pvtu in datasets
    )
    text = f"""<?xml version="1.0"?>
<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">
  <Collection>
{rows}
  </Collection>
</VTKFile>
"""
    path.write_text(text, encoding="utf-8")


def safe_time_label(index: int, value: str) -> str:
    label = re.sub(r"[^0-9A-Za-z_.+-]+", "_", value).strip("_")
    return f"t{index:04d}_{label}" if label else f"t{index:04d}"


def split_evenly(items: tuple[Piece, ...], group_count: int) -> list[tuple[Piece, ...]]:
    group_count = max(1, min(group_count, len(items)))
    base_size, remainder = divmod(len(items), group_count)
    groups: list[tuple[Piece, ...]] = []
    start = 0
    for group_index in range(group_count):
        end = start + base_size + (1 if group_index < remainder else 0)
        groups.append(items[start:end])
        start = end
    return groups


def validate_piece(piece: Piece) -> tuple[int, int]:
    dimensions = piece.topology_dimensions
    nodes = point_count(dimensions)
    if piece.geometry.dimensions != (nodes, 3):
        raise ValueError(
            f"grid {piece.name} geometry dimensions {piece.geometry.dimensions} "
            f"do not match topology {dimensions}"
        )
    if data_item_count(piece.attribute.dimensions) != nodes * piece.attribute_components:
        raise ValueError(
            f"grid {piece.name} attribute dimensions {piece.attribute.dimensions} "
            f"do not match topology {dimensions}"
        )
    if piece.geometry.precision != piece.attribute.precision:
        raise ValueError(f"grid {piece.name} mixes coordinate and value precision")
    return nodes, piece.attribute.precision


def write_piece_group_vtu(
    path: Path,
    pieces: tuple[Piece, ...],
    reader: BinaryItemReader,
    attribute_name: str,
    attribute_components: int,
) -> int:
    points = bytearray()
    values = bytearray()
    connectivity = bytearray()
    offsets = bytearray()
    types = bytearray()
    total_points = 0
    total_cells = 0
    precision = pieces[0].attribute.precision
    cell_offset = 0

    for piece in pieces:
        nodes, piece_precision = validate_piece(piece)
        if piece_precision != precision:
            raise ValueError("all pieces in a VTU file must use the same precision")
        if piece.attribute_components != attribute_components:
            raise ValueError("all pieces must use the same number of attribute components")

        piece_connectivity, piece_offsets, piece_types, piece_cells = make_cells(
            piece.topology_dimensions, total_points
        )
        points.extend(reader.read(piece.geometry))
        values.extend(reader.read(piece.attribute))
        connectivity.extend(piece_connectivity)
        for offset_index in range(piece_cells):
            (local_offset,) = struct.unpack_from("<i", piece_offsets, offset_index * 4)
            offsets.extend(struct.pack("<i", local_offset + cell_offset))
        types.extend(piece_types)

        if piece_offsets:
            (last_local_offset,) = struct.unpack_from("<i", piece_offsets, (piece_cells - 1) * 4)
            cell_offset += last_local_offset
        total_points += nodes
        total_cells += piece_cells

    write_vtu(
        path,
        points=bytes(points),
        values=bytes(values),
        cells=(bytes(connectivity), bytes(offsets), bytes(types)),
        point_count=total_points,
        cell_count=total_cells,
        attribute_name=attribute_name,
        attribute_components=attribute_components,
        precision=precision,
    )
    return precision


def convert(
    sdmf_file_path: Path,
    output_dir: Path,
    attribute_name: str | None,
    vtu_files_per_pvtu: int = DEFAULT_VTU_FILES_PER_PVTU,
) -> None:
    if vtu_files_per_pvtu < 1:
        raise ValueError("--vtu-files-per-pvtu must be at least 1")

    attribute_name, timesteps = parse_sdmf_file(sdmf_file_path, attribute_name)
    base_dir = sdmf_file_path.parent
    stem = sdmf_file_path.name
    if stem.endswith(".txt"):
        stem = stem[:-4]

    pvd_entries: list[tuple[str, Path]] = []

    with BinaryItemReader(base_dir) as reader:
        for timestep in timesteps:
            time_label = safe_time_label(timestep.index, timestep.value)
            timestep_dir = output_dir / time_label
            piece_paths: list[Path] = []
            timestep_precision = 4
            timestep_components = timestep.pieces[0].attribute_components

            for group_index, pieces in enumerate(split_evenly(timestep.pieces, vtu_files_per_pvtu)):
                piece_path = timestep_dir / f"{stem}_{time_label}_piece{group_index:04d}.vtu"
                timestep_precision = write_piece_group_vtu(
                    piece_path,
                    pieces,
                    reader,
                    attribute_name,
                    timestep_components,
                )
                piece_paths.append(piece_path)

            pvtu_path = output_dir / f"{stem}_{time_label}.pvtu"
            write_pvtu(pvtu_path, piece_paths, attribute_name, timestep_components, timestep_precision)
            pvd_entries.append((timestep.value, pvtu_path))

    write_pvd(output_dir / f"{stem}.pvd", pvd_entries)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert an SDMF 3DSMesh binary dataset to VTK PVTU/VTU files."
    )
    parser.add_argument("sdmf_file", type=Path, help="input SDMF file (has ending .txt)")
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        help="output directory; defaults to <input-stem>_pvtu beside the SDMF text file",
    )
    parser.add_argument(
        "-a",
        "--attribute",
        help="node attribute name to convert; defaults to the first node attribute",
    )
    parser.add_argument(
        "-n",
        "--vtu-files-per-pvtu",
        type=int,
        default=DEFAULT_VTU_FILES_PER_PVTU,
        help=f"maximum number of .vtu files referenced by each .pvtu (default: {DEFAULT_VTU_FILES_PER_PVTU})",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    sdmf_file_path = args.sdmf_file.resolve()
    output_dir = args.output_dir
    if output_dir is None:
        output_dir = sdmf_file_path.with_suffix("").parent / f"{sdmf_file_path.with_suffix('').name}_pvtu"
    output_dir = output_dir.resolve()

    try:
        convert(sdmf_file_path, output_dir, args.attribute, args.vtu_files_per_pvtu)
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    print(f"wrote {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

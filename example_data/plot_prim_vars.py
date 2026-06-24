#!/usr/bin/env python3
"""
Plot Parthenon primitive variables exposed by the openPMD VisIt plugin.

Run with:
  visit -cli -nowin -s example_data/plot_prim_vars.py /path/to/file.pmd

The script writes a metadata dump and one PNG per plotted variable.  By default
it plots every scalar whose metadata name starts with "prim/".
"""

import argparse
import os
import re
import sys


DEFAULT_PATTERNS = (r"^prim/",)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot prim/* variables from an openPMD .pmd file with VisIt."
    )
    parser.add_argument("dataset", help="Path to the .pmd descriptor to open.")
    parser.add_argument(
        "--pattern",
        action="append",
        default=None,
        help=(
            "Regex matched against VisIt scalar names. May be repeated. "
            "Default: '^prim/'."
        ),
    )
    parser.add_argument(
        "--var",
        action="append",
        default=None,
        help="Exact scalar variable name to plot. May be repeated.",
    )
    parser.add_argument(
        "--out-dir",
        default=None,
        help="Output directory for metadata and PNG files.",
    )
    parser.add_argument(
        "--format",
        default="PNG",
        choices=("PNG", "JPEG", "TIFF"),
        help="SaveWindow output format.",
    )
    parser.add_argument("--width", type=int, default=1400)
    parser.add_argument("--height", type=int, default=1000)
    parser.add_argument(
        "--list-only",
        action="store_true",
        help="Only print/dump metadata; do not create plots.",
    )
    return parser.parse_args()


def abs_existing_path(path):
    resolved = os.path.abspath(path)
    if not os.path.exists(resolved):
        raise RuntimeError("Dataset does not exist: {}".format(resolved))
    return resolved


def safe_name(name):
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", name.strip("/"))
    return cleaned or "unnamed"


def scalar_names(metadata):
    return [metadata.GetScalars(i).name for i in range(metadata.GetNumScalars())]


def vector_names(metadata):
    return [metadata.GetVectors(i).name for i in range(metadata.GetNumVectors())]


def mesh_names(metadata):
    return [metadata.GetMeshes(i).name for i in range(metadata.GetNumMeshes())]


def write_metadata_dump(path, metadata, scalars, vectors, meshes):
    with open(path, "w") as output:
        output.write("Meshes ({})\n".format(len(meshes)))
        for name in meshes:
            output.write("  {}\n".format(name))
        output.write("\nScalars ({})\n".format(len(scalars)))
        for name in scalars:
            output.write("  {}\n".format(name))
        output.write("\nVectors ({})\n".format(len(vectors)))
        for name in vectors:
            output.write("  {}\n".format(name))


def selected_scalars(args, scalars):
    selected = []

    if args.var:
        scalar_set = set(scalars)
        for name in args.var:
            if name in scalar_set:
                selected.append(name)
            else:
                print("Requested variable is not in VisIt scalar metadata: {}".format(name))
        return selected

    patterns = args.pattern if args.pattern else DEFAULT_PATTERNS
    compiled = [re.compile(pattern) for pattern in patterns]
    for name in scalars:
        if any(pattern.search(name) for pattern in compiled):
            selected.append(name)
    return selected


def configure_save_window(out_dir, file_stem, image_format, width, height):
    attrs = SaveWindowAttributes()
    attrs.outputToCurrentDirectory = 0
    attrs.outputDirectory = out_dir
    attrs.fileName = file_stem
    attrs.family = 0
    attrs.width = width
    attrs.height = height
    attrs.screenCapture = 0
    if image_format == "PNG":
        attrs.format = attrs.PNG
    elif image_format == "JPEG":
        attrs.format = attrs.JPEG
    elif image_format == "TIFF":
        attrs.format = attrs.TIFF
    SetSaveWindowAttributes(attrs)


def plot_scalar(var_name, out_dir, image_format, width, height):
    DeleteAllPlots()
    print("Plotting {}".format(var_name))
    AddPlot("Pseudocolor", var_name)
    DrawPlots()

    try:
        Query("MinMax")
        print("  MinMax: {}".format(GetQueryOutputString().strip()))
    except Exception as exc:
        print("  MinMax query failed: {}".format(exc))

    configure_save_window(out_dir, safe_name(var_name), image_format, width, height)
    saved = SaveWindow()
    print("  SaveWindow returned {}".format(saved))
    DeleteAllPlots()


def main():
    args = parse_args()
    dataset = abs_existing_path(args.dataset)
    out_dir = args.out_dir
    if out_dir is None:
        out_dir = os.path.join(os.getcwd(), "visit_prim_plots")
    out_dir = os.path.abspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)

    print("Opening dataset: {}".format(dataset))
    OpenDatabase(dataset)
    metadata = GetMetaData(dataset)

    meshes = mesh_names(metadata)
    scalars = scalar_names(metadata)
    vectors = vector_names(metadata)
    dump_path = os.path.join(out_dir, "metadata.txt")
    write_metadata_dump(dump_path, metadata, scalars, vectors, meshes)

    print("Metadata dump: {}".format(dump_path))
    print("Meshes: {}".format(len(meshes)))
    print("Scalars: {}".format(len(scalars)))
    print("Vectors: {}".format(len(vectors)))

    candidates = selected_scalars(args, scalars)
    print("Selected scalar vars: {}".format(len(candidates)))
    for name in candidates:
        print("  {}".format(name))

    if args.list_only:
        return 0

    if not candidates:
        print("No matching scalar variables were visible to VisIt.")
        print("Inspect {} to see the exact metadata names.".format(dump_path))
        return 2

    for name in candidates:
        try:
            plot_scalar(name, out_dir, args.format, args.width, args.height)
        except Exception as exc:
            print("FAILED to plot {}: {}".format(name, exc))

    CloseDatabase(dataset)
    return 0


if __name__ == "__main__":
    sys.exit(main())

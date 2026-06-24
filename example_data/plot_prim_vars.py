#!/usr/bin/env python3
"""
Plot Parthenon primitive variables exposed by the openPMD VisIt plugin.

Run with:
  visit -cli -nowin -s example_data/plot_prim_vars.py /path/to/file.pmd

The script writes a metadata dump and one PNG per plotted variable.  By default
it plots every scalar whose metadata name starts with "prim_".
"""

import argparse
import os
import re
import sys


DEFAULT_PATTERNS = (r"^prim_",)


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
            "Default: '^prim_'."
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
    parser.add_argument(
        "--dump-states",
        action="store_true",
        help="Dump metadata after walking every active time-slider state.",
    )
    parser.add_argument(
        "--plot-states",
        action="store_true",
        help="With --dump-states, plot selected variables at each state.",
    )
    parser.add_argument(
        "--state-limit",
        type=int,
        default=None,
        help="Maximum number of time-slider states to dump in --dump-states mode.",
    )
    parser.add_argument(
        "--initial-state",
        type=int,
        default=None,
        help="Set the active time-slider state before listing or plotting.",
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
        for i, name in enumerate(meshes):
            mesh = metadata.GetMeshes(i)
            output.write("  {}\n".format(name))
            for attr in (
                "meshType",
                "topologicalDimension",
                "spatialDimension",
                "numBlocks",
                "numGroups",
                "validVariable",
                "hideFromGUI",
            ):
                if hasattr(mesh, attr):
                    output.write("    {}={}\n".format(attr, getattr(mesh, attr)))
        output.write("\nScalars ({})\n".format(len(scalars)))
        for i, name in enumerate(scalars):
            scalar = metadata.GetScalars(i)
            output.write("  {}\n".format(name))
            for attr in (
                "originalName",
                "meshName",
                "centering",
                "validVariable",
                "hideFromGUI",
                "treatAsASCII",
            ):
                if hasattr(scalar, attr):
                    output.write("    {}={}\n".format(attr, getattr(scalar, attr)))
        output.write("\nVectors ({})\n".format(len(vectors)))
        for i, name in enumerate(vectors):
            vector = metadata.GetVectors(i)
            output.write("  {}\n".format(name))
            for attr in (
                "originalName",
                "meshName",
                "centering",
                "validVariable",
                "hideFromGUI",
                "varDim",
            ):
                if hasattr(vector, attr):
                    output.write("    {}={}\n".format(attr, getattr(vector, attr)))


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


def dump_metadata_for_state(args, dataset, state, out_dir):
    print("Reading metadata for time state {}".format(state))
    if args.plot_states:
        result = TimeSliderSetState(state)
        print("  TimeSliderSetState returned {}".format(result))
        metadata = GetMetaData(dataset)
    else:
        metadata = GetMetaData(dataset, state)
    meshes = mesh_names(metadata)
    scalars = scalar_names(metadata)
    vectors = vector_names(metadata)
    dump_path = os.path.join(out_dir, "metadata_state_{:05d}.txt".format(state))
    write_metadata_dump(dump_path, metadata, scalars, vectors, meshes)
    candidates = selected_scalars(args, scalars)
    print(
        "  state {}: meshes={} scalars={} vectors={} prim_scalars={}".format(
            state, len(meshes), len(scalars), len(vectors), len(candidates)
        )
    )
    print("  metadata dump: {}".format(dump_path))
    for name in candidates:
        print("    {}".format(name))
    if args.plot_states:
        state_plot_dir = os.path.join(out_dir, "state_{:05d}".format(state))
        os.makedirs(state_plot_dir, exist_ok=True)
        for name in candidates:
            plot_scalar(
                name, state_plot_dir, args.format, args.width, args.height
            )
    return candidates


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
    if args.initial_state is not None:
        print("Setting initial time state: {}".format(args.initial_state))
        result = TimeSliderSetState(args.initial_state)
        print("TimeSliderSetState returned {}".format(result))

    if args.dump_states:
        num_states = TimeSliderGetNStates()
        print("Time slider states: {}".format(num_states))
        missing_states = []
        state_count = num_states
        if args.state_limit is not None:
            state_count = min(num_states, max(0, args.state_limit))
        for state in range(state_count):
            candidates = dump_metadata_for_state(args, dataset, state, out_dir)
            if not candidates:
                missing_states.append(state)
        if missing_states:
            print("No prim scalars were visible in states: {}".format(missing_states))
            return 2
        return 0

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

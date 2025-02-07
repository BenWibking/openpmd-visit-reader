# [WIP] openpmd-visit-reader

This is a serial database plugin for VisIt that reads files using the OpenPMD-api.

Supported:
* Cartesian geometry
* Uniform grids
* Cell-centered variables
* Node-centered variables
* Particles

Buggy:
* Data layout transformation
  * Temporary workaround:
    * Override mesh axis labels (specify new labels on line 2 of *.pmd file): https://github.com/BenWibking/openpmd-visit-reader/blob/3c8e5dcfbfb1169ab21d6ee9c47970ee6a5e71a5/example_data/hdf5_2d.pmd#L2
    * Override particle axis labels (specify new labels on line 3 of *.pmd file): https://github.com/BenWibking/openpmd-visit-reader/blob/3c8e5dcfbfb1169ab21d6ee9c47970ee6a5e71a5/example_data/hdf5_2d.pmd#L3

Currently unsupported:
* Particle data (other than positions)
* AMR meshes
  * Used in [WarpX](https://github.com/ECP-WarpX/WarpX/blob/d79fe71ae810364b02017ef70c82c70f667c8e19/Source/Diagnostics/WarpXOpenPMD.cpp#L1282). Each AMR level is named `*_lvlN`. Each AMR patch is saved as a [separate ADIOS2 chunk](https://github.com/ECP-WarpX/WarpX/blob/d79fe71ae810364b02017ef70c82c70f667c8e19/Source/Diagnostics/WarpXOpenPMD.cpp#L1462). Reconstructing the AMR mesh is up to the reader.
* Parallel reads using MPI

Unplanned:
* non-Cartesian geometry
* Face-centered variables
* Edge-centered variables

## Building

1. First, make sure all of the submodules are up to date:
   ```
   git submodule update --init
   ```

2. Edit this line to point to your home directory:
https://github.com/BenWibking/openpmd-visit-reader/blob/856c967442b6b903f3ddac2fc49963b40f547169/CMakeLists.txt#L3

3. Edit this line to point to your VisIt installation:
https://github.com/BenWibking/openpmd-visit-reader/blob/856c967442b6b903f3ddac2fc49963b40f547169/CMakeLists.txt#L4

4. Copy `VisItLibraryDependencies.cmake` from your VisIt installation to the root of the repository.

5. *(macOS only)* Run:
   ```
   ./fix_VisItLibraryDependencies_macos.sh
   ```

6. Finally, build with:
   ```
   mkdir build && cd build
   cmake .. -GNinja
   ninja
   ```
   The plugin should be installed to your `~/.visit` directory, where VisIt should detect and load it automatically. After recompiling the plugin, you may have to restart VisIt in order to use the new version of the plugin.

## Testing

1. Extract the example data using the provided script:
   ```
   cd example_data
   ./extract_example_data.sh
   ```

2. Load any of the `*.pmd` files in the `example_data/` directory in VisIt (clicking `Open` in the main window). This file extension is associated with OpenPMD files, following the [convention](https://openpmd-api.readthedocs.io/en/latest/analysis/paraview.html#openpmd) used for the Paraview reader plugin.

   As an extension specific to this plugin, lines 2 and 3 of the `*.pmd` file can be used to override the as-written axis labels for mesh and particle data. This is used to avoid data layout transforms, which do not currently work correctly.

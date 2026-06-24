# openpmd-visit-reader (beta)

This is a parallel database plugin for VisIt that reads files using the OpenPMD-api.

Supported:
* Cartesian geometry
* AMR meshes
  * Used in [WarpX](https://github.com/ECP-WarpX/WarpX/blob/d79fe71ae810364b02017ef70c82c70f667c8e19/Source/Diagnostics/WarpXOpenPMD.cpp#L1282). Each AMR level is named `*_lvlN`. Each AMR patch is saved as a [separate ADIOS2 chunk](https://github.com/ECP-WarpX/WarpX/blob/d79fe71ae810364b02017ef70c82c70f667c8e19/Source/Diagnostics/WarpXOpenPMD.cpp#L1462) with level-centric patch extents.
* Cell-centered scalar variables
* Cell-centered vector variables
* Particle species (rendered as point meshes with per-particle scalars and vectors)

Buggy:
* Node-centered variables (see https://github.com/BenWibking/openpmd-visit-reader/issues/1)

Unplanned:
* non-Cartesian geometry
* Face-centered variables
* Edge-centered variables

## Building

1. First, make sure all of the submodules are up to date:
   ```
   git submodule update --init
   ```

2. Build with (edit `VISIT_PLUGIN_DIR` and `VISIT_PLUGIN_VS_INSTALL_FILE` according to your installation):
   ```
   mkdir build && cd build
   cmake .. -GNinja -DVISIT_PLUGIN_DIR="${HOME}/.visit/3.5.0/darwin-arm64/plugins" -DVISIT_PLUGIN_VS_INSTALL_FILE="/Applications/VisIt.app/Contents/Resources/3.5.0/darwin-arm64/include/PluginVsInstall.cmake"
   ninja
   ```
   The plugin should be installed to your `~/.visit` directory, where VisIt should detect and load it automatically. After recompiling the plugin, you may have to restart VisIt in order to use the new version of the plugin.

   `CMakeLists.txt` automatically reads `VisItLibraryDependencies.cmake` from the same VisIt include directory as `PluginVsInstall.cmake`. On macOS, it also normalizes VisIt's generated VTK dependency names to the shipped dylib filenames, so no manual dependency-file patching is needed.

   Keep the `IopenpmdDatabase` info plugin lightweight. VisIt loads every database info plugin while populating the `Open` dialog, before any OpenPMD file is selected. If `libIopenpmdDatabase.dylib` links against OpenPMD, ADIOS2, or other backend dependencies, a loader/signing/rpath failure in those libraries can prevent `mdserver` from returning the directory listing and leave the `Open` dialog empty. The OpenPMD dependencies should stay in the metadata-server and engine plugin libraries that are used when a `.pmd` file is actually opened.

   On macOS, you can also use:
   ```
   ./build_macos.sh
   ```

   To skip caching the structured domain boundary metadata (and therefore disable VisIt's ghost-synthesis metadata path), add `-DOPENPMD_DISABLE_STRUCTURED_BOUNDARY_CACHE=ON` when invoking CMake.

## Testing

1. Extract the example data using the provided script:
   ```
   cd example_data
   ./extract_example_data.sh
   ```

2. Load any of the `*.pmd` files in the `example_data/` directory in VisIt (clicking `Open` in the main window). This file extension is associated with OpenPMD files, following the [convention](https://openpmd-api.readthedocs.io/en/latest/analysis/paraview.html#openpmd) used for the Paraview reader plugin.
3. For datasets that define particle species, select the `<species>_particles` mesh in the `Meshes` pane and plot associated scalars/vectors to confirm coordinate loading and per-particle data lookups.

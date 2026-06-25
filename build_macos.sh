## do a standard build with VisIt's bundled MPI
set -x

cmake -S . -B build -GNinja -DopenPMD_USE_ADIOS2=ON -DopenPMD_USE_HDF5=ON -DopenPMD_USE_MPI=ON -DVISIT_PLUGIN_DIR="${HOME}/.visit/3.5.0/darwin-arm64/plugins" -DVISIT_PLUGIN_VS_INSTALL_FILE="/Applications/VisIt.app/Contents/Resources/3.5.0/darwin-arm64/include/PluginVsInstall.cmake" "$@"

cmake --build build
# the plugin should actually be installed to ~/.visit in the build step
#cmake --install build

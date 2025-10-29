## do a standard build without MPI
set -x

cmake -S . -B build -GNinja -DopenPMD_USE_ADIOS2=ON -DopenPMD_USE_HDF5=ON -DopenPMD_USE_MPI=OFF -DVISIT_PLUGIN_DIR="${HOME}/.visit/3.4.2/darwin-arm64/plugins" -DVISIT_PLUGIN_VS_INSTALL_FILE="/Applications/VisIt.app/Contents/Resources/3.4.2/darwin-arm64/include/PluginVsInstall.cmake"

cmake --build build

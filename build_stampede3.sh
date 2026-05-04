#!/usr/bin/env bash
set -euo pipefail
set -x

cmake -S . -B build -GNinja \
  -DVISIT_PLUGIN_DIR="${HOME}/.visit/3.5.0/linux-x86_64/plugins" \
  -DVISIT_PLUGIN_VS_INSTALL_FILE="${HOME}/visit3_5_0.linux-x86_64/3.5.0/linux-x86_64/include/PluginVsInstall.cmake" \
  -DopenPMD_USE_ADIOS2=ON \
  -DopenPMD_USE_HDF5=OFF \
  -DopenPMD_USE_MPI=OFF \
  -DOPENPMD_USE_VISIT_ADIOS2=ON \
  "$@"

cmake --build build -j16

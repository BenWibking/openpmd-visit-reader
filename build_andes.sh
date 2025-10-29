## build with ADIOS2 + MPI support on Andes
## NOTE: you can't use any non-default modules on Andes -- they don't get loaded when Visit is run in client-server mode
module reset

# build Blosc2
# ...

# build ADIOS2 with Blosc2 support
# ...

set -x
#CXX=g++ CC=gcc ## DO NOT SET
cmake -S . -B build -DopenPMD_USE_ADIOS2=ON -DopenPMD_USE_HDF5=OFF -DopenPMD_USE_MPI=ON
cmake --build build -j16

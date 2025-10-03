## do a standard build without MPI
set -x
cmake -S . -B build -GNinja -DopenPMD_USE_ADIOS2=ON -DopenPMD_USE_HDF5=ON -DopenPMD_USE_MPI=OFF
cmake --build build

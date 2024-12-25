# openpmd-visit-reader

This is a serial database plugin for VisIt that reads files using the OpenPMD-api.

Supported:
* Cartesian geometry
* Uniform grids
* Cell-centered variables
* Node-centered variables
* Particles

Buggy:
* Data layout transformation

Currently unsupported:
* Particle data (other than positions)
* AMR meshes
  * Used in [WarpX](https://github.com/ECP-WarpX/WarpX/blob/d79fe71ae810364b02017ef70c82c70f667c8e19/Source/Diagnostics/WarpXOpenPMD.cpp#L1282). Each AMR level is named `*_lvlN`. Each AMR patch is saved as a [separate ADIOS2 chunk](https://github.com/ECP-WarpX/WarpX/blob/d79fe71ae810364b02017ef70c82c70f667c8e19/Source/Diagnostics/WarpXOpenPMD.cpp#L1462). Reconstructing the AMR mesh is up to the reader.
* Parallel reads using MPI

Unplanned:
* non-Cartesian geometry
* Face-centered variables
* Edge-centered variables

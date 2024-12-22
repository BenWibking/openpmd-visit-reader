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
* AMR meshes
* Parallel reads using MPI

Unplanned:
* non-Cartesian geometry
* Face-centered variables
* Edge-centered variables

# openpmd-visit-reader

This is a serial database plugin for VisIt that reads files using the OpenPMD-api.

Supported:
* Cartesian geometry
* Uniform grids
* Cell-centered variables

Not working:
* Node-centered variables

Currently unsupported:
* Particles
* AMR meshes
* non-Cartesian geometry
* Face-, edge-centered variables
* Parallel reads (using MPI)

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

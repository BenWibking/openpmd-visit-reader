// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ****************************************************************************
//  avtopenpmdFileFormat.h
// ****************************************************************************

#ifndef AVT_openpmd_FILE_FORMAT_H
#define AVT_openpmd_FILE_FORMAT_H

#include <string>
#include <unordered_map>

#include <openPMD/openPMD.hpp>

#include <DebugStream.h>
#include <avtMTSDFileFormat.h>

struct GeometryData {
  std::vector<std::string> axisLabels;
  std::vector<uint64_t> extent;
  std::vector<double> gridSpacing;
  std::vector<double> gridOrigin;
};

// ****************************************************************************
//  Class: avtopenpmdFileFormat
//
//  Purpose:
//      Reads in openpmd files as a plugin to VisIt.
//
//  Programmer: benwibking -- generated by xml2avt
//  Creation:   Fri Dec 6 17:16:49 PST 2024
//
// ****************************************************************************

class avtopenpmdFileFormat : public avtMTSDFileFormat {
public:
  avtopenpmdFileFormat(const char *);
  virtual ~avtopenpmdFileFormat() { ; }

  //
  // This is used to return unconvention data -- ranging from material
  // information to information about block connectivity.
  //
  // virtual void      *GetAuxiliaryData(const char *var, int timestep,
  //                                     const char *type, void *args,
  //                                     DestructorFunction &);
  //

  //
  // If you know the times and cycle numbers, overload this function.
  // Otherwise, VisIt will make up some reasonable ones for you.
  //
  // virtual void        GetCycles(std::vector<int> &);
  // virtual void        GetTimes(std::vector<double> &);
  //

  virtual int GetNTimesteps(void);

  virtual const char *GetType(void) { return "openpmd"; }
  virtual void FreeUpResources(void);

  virtual vtkDataSet *GetMesh(int, const char *);
  virtual vtkDataArray *GetVar(int, const char *);
  virtual vtkDataArray *GetVectorVar(int, const char *);

  enum class DatasetType { Field = 0, ParticleSpecies };

  // Order of transposes relative to the layout of data
  // For 3D data : in case of (z, y, x) ordering of axes
  const std::vector<int> _3d = {2, 1, 0};
  // For 2D data : in case of (z, x) ordering of axes
  const std::vector<int> _2d_xz = {1, 2, 0};
  // For 2D data : in case of (y, x) ordering of axes
  const std::vector<int> _2d_xy = {1, 0, 2};
  // For 1D data : in case data layout (z)
  const std::vector<int> _1d = {2, 1, 0};
  // identity transpose
  const std::vector<int> _identity = {0, 1, 2};

protected:
  // DATA MEMBERS
  openPMD::Series series_;
  std::vector<unsigned long long> iterationIndex_;
  std::unordered_map<std::string, std::tuple<std::string, std::string>>
      varMap_; // from VisIt varname, get openPMD mesh name AND record component
               // name
  std::unordered_map<std::string, std::tuple<DatasetType, std::string>>
      meshMap_; // from VisIt mesh name, get openPMD mesh name AND DatasetType
  bool doOverrideAxisOrder_{false};
  std::vector<std::string> overrideAxisLabels_;

  virtual void PopulateDatabaseMetaData(avtDatabaseMetaData *, int);
  void ReadFieldMetaData(avtDatabaseMetaData *md, openPMD::Iteration const &i);
  void ReadParticleMetaData(avtDatabaseMetaData *md,
                            openPMD::Iteration const &i);

  vtkDataSet *GetMeshField(openPMD::Iteration i, std::string const &meshname);
  vtkDataSet *GetMeshParticles(openPMD::Iteration i,
                               std::string const &meshname);

  template <typename T> void ScaleVarData(T *xyz_ptr, size_t nelem, T unitSI);

  template <typename T> avtCentering GetCenteringType(T const &mesh);

  GeometryData GetGeometry3D(openPMD::Mesh const &mesh);
  GeometryData GetGeometryXYZ(openPMD::Mesh const &mesh);

  std::vector<int> GetIndexOrder(std::vector<std::string> axisLabels);

  std::vector<int> GetAxisTranspose(std::vector<std::string> const &axisLabels);

  template <typename T>
  void TransposeVector(std::vector<T> &vec_to_transpose,
                       std::vector<int> const &transpose);

  std::tuple<size_t, size_t, size_t>
  GetIndexCoefficients(std::vector<int> const &indexOrder,
                       std::vector<uint64_t> const &ndims);

  template <typename T>
  void TransposeArray(T *data_ptr, openPMD::Mesh const &mesh,
                      openPMD::Mesh::MeshRecordComponent const &rcomp);
};

template <typename T>
avtCentering avtopenpmdFileFormat::GetCenteringType(T const &mesh) {
  // read the element centering
  std::vector<float> centering;
  try {
    centering = mesh.template position<float>();
  } catch (openPMD::Error) {
    debug5 << "[openpmd-api-plugin] "
           << "Can't read centering for mesh!\n";
    return AVT_UNKNOWN_CENT;
  }

  // cell-centered == {0.5, 0.5, 0.5}
  bool isCellCentered = true;
  for (int idim = 0; idim < centering.size(); idim++) {
    isCellCentered = (isCellCentered && (centering[idim] == 0.5));
  }

  // node-centered == {0, 0, 0}
  bool isNodeCentered = true;
  for (int idim = 0; idim < centering.size(); idim++) {
    isNodeCentered = (isNodeCentered && (centering[idim] == 0.));
  }

  if (isCellCentered) {
    debug5 << "[openpmd-api-plugin] "
           << "Mesh is cell-centered.\n";
    return AVT_ZONECENT;
  } else if (isNodeCentered) {
    debug5 << "[openpmd-api-plugin] "
           << "Mesh is node-centered.\n";
    return AVT_NODECENT;
  }

  // face-centered == {0, 0.5, 0.5}, etc.
  // edge-centered == {0, 0, 0.5}, etc.
  // However, all other centerings are unsupported by VisIt.
  debug5 << "[openpmd-api-plugin] "
         << "Mesh has unsupported centering.\n";
  return AVT_UNKNOWN_CENT;
}

template <typename T>
void avtopenpmdFileFormat::ScaleVarData(T *xyz_ptr, size_t nelem, T unitSI) {
  for (size_t idx = 0; idx < nelem; idx++) {
    xyz_ptr[idx] *= unitSI;
  }
}

#endif

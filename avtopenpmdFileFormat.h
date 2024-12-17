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

#include <avtMTSDFileFormat.h>

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

protected:
  // DATA MEMBERS
  openPMD::Series series_;
  std::vector<unsigned long long> iterationIndex_;
  std::unordered_map<std::string, std::tuple<std::string, std::string>>
      varMap_; // from VisIt varname, get openPMD mesh name AND record component
               // name
  std::unordered_map<std::string, std::tuple<DatasetType, std::string>>
      meshMap_; // from VisIt mesh name, get openPMD mesh name AND DatasetType

  virtual void PopulateDatabaseMetaData(avtDatabaseMetaData *, int);
  void ReadFieldMetaData(avtDatabaseMetaData *md, openPMD::Iteration const &i);
  void ReadParticleMetaData(avtDatabaseMetaData *md,
                            openPMD::Iteration const &i);

  vtkDataSet *GetMeshField(openPMD::Iteration i, std::string const &meshname);
  vtkDataSet *GetMeshParticles(openPMD::Iteration i,
                               std::string const &meshname);

  template <typename T> avtCentering GetCenteringType(T const &mesh);
};

#endif

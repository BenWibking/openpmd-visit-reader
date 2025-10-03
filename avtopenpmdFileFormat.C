// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ****************************************************************************
//  avtopenpmdFileFormat.C
// ****************************************************************************

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <csignal>
#include <atomic>
#include <iostream>
#include <map>
#include <limits>
#include <sstream>
#include <cstring>
#include <unordered_map>
#include <stdexcept>
#include <string>
#include <unistd.h>

#include <avtDatabaseMetaData.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkRectilinearGrid.h>

#include <openPMD/openPMD.hpp>

#include <cpptrace/cpptrace.hpp>

#include <DBOptionsAttributes.h>
#include <DebugStream.h>
#include <Expression.h>
#include <InvalidVariableException.h>
#include <InvalidFilesException.h>

#include "avtopenpmdFileFormat.h"
#include "dataLayoutTransform.h" // NOLINT(unused-includes)

namespace {

inline void ComputeLogicalExtents(PatchInfo &patch) {
  for (int axis = 0; axis < 3; ++axis) {
    int lower = 0;
    if (axis < static_cast<int>(patch.offset.size())) {
      lower = static_cast<int>(patch.offset[axis]);
    }

    int cells = 1;
    if (axis < static_cast<int>(patch.extent.size())) {
      cells = static_cast<int>(patch.extent[axis]);
    }
    if (patch.centering == AVT_NODECENT && cells > 1) {
      cells -= 1;
    }
    if (cells <= 0) {
      cells = 1;
    }

    patch.logicalLower[axis] = lower;
    patch.logicalUpper[axis] = lower + cells - 1;
  }
}

inline void GetPhysicalBounds(const PatchInfo &patch, double minBounds[3],
                              double maxBounds[3]) {
  for (int axis = 0; axis < 3; ++axis) {
    const double spacing = patch.spacing[axis];
    const double origin = patch.origin[axis];
    const int cells = patch.logicalUpper[axis] - patch.logicalLower[axis] + 1;
    if (spacing == 0.0 || cells <= 0) {
      minBounds[axis] = origin;
      maxBounds[axis] = origin;
    } else {
      minBounds[axis] = origin;
      maxBounds[axis] = origin + spacing * static_cast<double>(cells);
    }
  }
}

inline bool IsChildPatch(const PatchInfo &coarse, const PatchInfo &fine) {
  double coarseMin[3] = {0.0, 0.0, 0.0};
  double coarseMax[3] = {0.0, 0.0, 0.0};
  double fineMin[3] = {0.0, 0.0, 0.0};
  double fineMax[3] = {0.0, 0.0, 0.0};

  GetPhysicalBounds(coarse, coarseMin, coarseMax);
  GetPhysicalBounds(fine, fineMin, fineMax);

  for (int axis = 0; axis < 3; ++axis) {
    double tol = std::max({std::abs(coarse.spacing[axis]),
                           std::abs(fine.spacing[axis]), 1.0}) * 1e-5;
    if (fineMin[axis] < coarseMin[axis] - tol) {
      return false;
    }
    if (fineMax[axis] > coarseMax[axis] + tol) {
      return false;
    }
  }
  return true;
}

inline std::vector<int> MakeRefinementVector(const std::array<int, 3> &ratio,
                                             int dims) {
  std::vector<int> result(std::max(1, dims), 1);
  for (int axis = 0; axis < std::min(3, static_cast<int>(result.size()));
       ++axis) {
    result[axis] = ratio[axis];
  }
  return result;
}

inline std::vector<double> MakeCellSizeVector(const std::array<double, 3> &sizes,
                                              int dims) {
  std::vector<double> result(std::max(1, dims), 0.0);
  for (int axis = 0; axis < std::min(3, static_cast<int>(result.size()));
       ++axis) {
    result[axis] = sizes[axis];
  }
  return result;
}

inline void WriteStderr(const char *message, size_t length) {
  if (message == nullptr || length == 0) {
    return;
  }
  ssize_t remaining = static_cast<ssize_t>(length);
  const char *ptr = message;
  while (remaining > 0) {
    ssize_t written = ::write(STDERR_FILENO, ptr, static_cast<size_t>(remaining));
    if (written <= 0) {
      break;
    }
    remaining -= written;
    ptr += written;
  }
}

std::atomic<bool> &SegvHandlerActiveFlag() {
  static std::atomic<bool> flag{false};
  return flag;
}

void PrintSegvStackTrace() {
  constexpr char header[] =
      "\n[openpmd-api-plugin] cpptrace captured stack trace (signal handler):\n";
  WriteStderr(header, sizeof(header) - 1);
  try {
    auto trace = cpptrace::generate_trace();
    std::string traceString = trace.to_string();
    WriteStderr(traceString.c_str(), traceString.size());
    const char newline = '\n';
    WriteStderr(&newline, 1);
  } catch (...) {
    constexpr char failure[] =
        "[openpmd-api-plugin] cpptrace failed to generate stack trace after signal.\n";
    WriteStderr(failure, sizeof(failure) - 1);
  }
}

void SegvSignalHandler(int sig, siginfo_t *, void *) {
  auto &active = SegvHandlerActiveFlag();
  bool expected = false;
  if (active.compare_exchange_strong(expected, true)) {
    PrintSegvStackTrace();
  }

  struct sigaction defaultAction;
  std::memset(&defaultAction, 0, sizeof(defaultAction));
  defaultAction.sa_handler = SIG_DFL;
  sigemptyset(&defaultAction.sa_mask);
  sigaction(sig, &defaultAction, nullptr);
  raise(sig);
}

void InstallSegfaultTraceHandler() {
  static std::atomic<bool> installed{false};
  bool expected = false;
  if (!installed.compare_exchange_strong(expected, true)) {
    return;
  }

  cpptrace::register_terminate_handler();
  (void)cpptrace::generate_trace();

  struct sigaction action;
  std::memset(&action, 0, sizeof(action));
  action.sa_sigaction = &SegvSignalHandler;
  sigemptyset(&action.sa_mask);
  action.sa_flags = SA_SIGINFO;
  if (sigaction(SIGSEGV, &action, nullptr) != 0) {
    constexpr char failure[] =
        "[openpmd-api-plugin] Failed to register cpptrace SIGSEGV handler.\n";
    WriteStderr(failure, sizeof(failure) - 1);
  }
}

void DeleteStructuredDomainNesting(void *ptr) {
  delete static_cast<avtStructuredDomainNesting *>(ptr);
}

void DeleteStructuredDomainBoundaries(void *ptr) {
  delete static_cast<avtStructuredDomainBoundaries *>(ptr);
}

template <typename Container>
inline std::string JoinContainer(const Container &values) {
  std::ostringstream oss;
  oss << '[';
  for (size_t idx = 0; idx < values.size(); ++idx) {
    if (idx > 0) {
      oss << ',';
    }
    oss << values[idx];
  }
  oss << ']';
  return oss.str();
}

template <typename T, size_t N>
inline std::string JoinArray(const T (&values)[N]) {
  std::ostringstream oss;
  oss << '[';
  for (size_t idx = 0; idx < N; ++idx) {
    if (idx > 0) {
      oss << ',';
    }
    oss << values[idx];
  }
  oss << ']';
  return oss.str();
}

inline const char *CenteringToString(avtCentering cent) {
  switch (cent) {
  case AVT_NODECENT:
    return "node";
  case AVT_ZONECENT:
    return "zone";
  case AVT_UNKNOWN_CENT:
    return "unknown";
#ifdef AVT_NO_CENTERING
  case AVT_NO_CENTERING:
    return "none";
#endif
  default:
    return "other";
  }
}

inline const char *DataOrderToString(openPMD::Mesh::DataOrder order) {
  switch (order) {
  case openPMD::Mesh::DataOrder::C:
    return "C";
  case openPMD::Mesh::DataOrder::F:
    return "F";
  default:
    return "unknown";
  }
}

inline void LogPatchSummary(const PatchInfo &patch,
                            const std::string &context) {
  debug3 << "[openpmd-api-plugin] " << context
         << " mesh=" << patch.meshName << " level=" << patch.level
         << " centering=" << CenteringToString(patch.centering)
         << " offset=" << JoinContainer(patch.offset)
         << " extent=" << JoinContainer(patch.extent)
         << " origin=" << JoinArray(patch.origin)
         << " spacing=" << JoinArray(patch.spacing)
         << " logicalLower=" << JoinArray(patch.logicalLower)
         << " logicalUpper=" << JoinArray(patch.logicalUpper) << '\n';
}

inline std::string JoinStrings(const std::vector<std::string> &values) {
  std::ostringstream oss;
  oss << '[';
  for (size_t idx = 0; idx < values.size(); ++idx) {
    if (idx > 0) {
      oss << ", ";
    }
    oss << values[idx];
  }
  oss << ']';
  return oss.str();
}

} // namespace

// ****************************************************************************
//  Method: avtopenpmdFileFormat constructor
//
//  Programmer: benwibking -- generated by xml2avt
//  Creation:   Fri Dec 6 17:16:49 PST 2024
//
// ****************************************************************************

avtopenpmdFileFormat::avtopenpmdFileFormat(const char *filename)
    : avtMTMDFileFormat(filename) {
  InstallSegfaultTraceHandler();
  debug1 << "[openpmd-api-plugin] Constructing reader for descriptor '"
         << filename << "'\n";
  //
  // initialize an openPMD::series object
  //

  // read incomplete filepath string from file 'filename'
  std::string opmd_filestring;
  std::string opmd_overrideMeshAxisLabels;
  std::string opmd_overrideParticleAxisLabels;

  {
    std::ifstream file(filename);
    if (!file.is_open()) {
      debug1 << "[openpmd-api-plugin] Failed to open descriptor file: "
             << filename << "\n";
    }

    // get file string
    std::getline(file, opmd_filestring);

    // if it exists, get overrideMeshAxisLabels
    std::getline(file, opmd_overrideMeshAxisLabels);
    if (opmd_overrideMeshAxisLabels != "") {
      doOverrideMeshAxisOrder_ = true;
      auto iss = std::istringstream{opmd_overrideMeshAxisLabels};
      auto str = std::string{};
      while (iss >> str) {
        overrideMeshAxisLabels_.push_back(str);
      }
      debug3 << "[openpmd-api-plugin] overrideMeshAxisLabels enabled: "
             << JoinStrings(overrideMeshAxisLabels_) << "\n";
    } else {
      debug3 << "[openpmd-api-plugin] overrideMeshAxisLabels disabled\n";
    }

    // if it exists, get overrideParticleAxisLabels
    std::getline(file, opmd_overrideParticleAxisLabels);
    if (opmd_overrideParticleAxisLabels != "") {
      doOverrideParticleAxisOrder_ = true;
      auto iss = std::istringstream{opmd_overrideParticleAxisLabels};
      auto str = std::string{};
      while (iss >> str) {
        overrideParticleAxisLabels_.push_back(str);
      }
      debug3 << "[openpmd-api-plugin] overrideParticleAxisLabels enabled: "
             << JoinStrings(overrideParticleAxisLabels_) << "\n";
    } else {
      debug3 << "[openpmd-api-plugin] overrideParticleAxisLabels disabled\n";
    }

    // close file
    file.close();
  }

  // construct complete filepath
  std::filesystem::path p(filename);
  std::string os_pathsep = "/";
  std::string parent_path = p.parent_path();
  std::string opmd_filepath = parent_path + os_pathsep + opmd_filestring;
  debug1 << "[openpmd-api-plugin] Resolved openPMD path: " << opmd_filepath
         << "\n";
  debug5 << "[openpmd-api-plugin] "
         << "Reading OpenPMD series: " << opmd_filepath << "\n";

  // open openPMD series
  try {
    series_ = openPMD::Series(opmd_filepath, openPMD::Access::READ_ONLY);
  } catch (std::exception const &ex) {
    std::string msg = std::string("Failed to open openPMD series '") +
                      opmd_filepath + "': " + ex.what();
    debug1 << "[openpmd-api-plugin] " << msg << "\n";
    EXCEPTION1(InvalidFilesException, msg.c_str());
  }
  debug5 << "[openpmd-api-plugin] "
         << "This file uses openPMD-standard version " << series_.openPMD()
         << '\n';

  debug1 << "[openpmd-api-plugin] Series openPMD version " << series_.openPMD()
         << " iterations=" << series_.snapshots().size() << "\n";

  // read iteration count
  debug5 << "[openpmd-api-plugin] "
         << "This file contains " << series_.snapshots().size()
         << " iterations:";
  const size_t numTimesteps = series_.snapshots().size();
  iterationIndex_ = std::vector<unsigned long long>(numTimesteps);
  meshHierarchyCache_.resize(numTimesteps);

  // save map from timeState to iteration index
  // NOTE: openPMD's iteration index can be an arbitrary number, and can
  // skip numbers. For instance, a dataset can have iterations = {550, 600}.
  int timeState = 0;
  for (auto const &iter : series_.snapshots()) {
    debug5 << "\n\t" << iter.first;
    iterationIndex_.at(timeState) = iter.first;
    timeState++;
  }
  debug5 << '\n';
  debug1 << "[openpmd-api-plugin] Iteration indices loaded for "
         << iterationIndex_.size() << " timesteps\n";

  debug1 << "[openpmd-api-plugin] Constructor about to finish. this="
         << static_cast<const void *>(this)
         << " cacheSize=" << meshHierarchyCache_.size()
         << " iterationsAvailable=" << series_.snapshots().size() << "\n";
}

avtopenpmdFileFormat::~avtopenpmdFileFormat() {
  debug1 << "[openpmd-api-plugin] Destructor invoked. this="
         << static_cast<const void *>(this) << "\n";
}

// ****************************************************************************
//  Method: avtopenpmdFileFormat::GetNTimesteps
//
//  Purpose:
//      Tells the rest of the code how many timesteps there are in this file.
//
//  Programmer: benwibking -- generated by xml2avt
//  Creation:   Fri Dec 6 17:16:49 PST 2024
//
// ****************************************************************************

int avtopenpmdFileFormat::GetNTimesteps(void) {
  return series_.snapshots().size();
}

// ****************************************************************************
//  Method: avtopenpmdFileFormat::FreeUpResources
//
//  Purpose:
//      When VisIt is done focusing on a particular timestep, it asks that
//      timestep to free up any resources (memory, file descriptors) that
//      it has associated with it.  This method is the mechanism for doing
//      that.
//
//  Programmer: benwibking -- generated by xml2avt
//  Creation:   Fri Dec 6 17:16:49 PST 2024
//
// ****************************************************************************

void avtopenpmdFileFormat::FreeUpResources(void) {
  debug1 << "[openpmd-api-plugin] FreeUpResources\n";
}

void avtopenpmdFileFormat::BuildFieldHierarchy(avtDatabaseMetaData *md,
                                               openPMD::Iteration const &i,
                                               int timeState) {
  debug1 << "[openpmd-api-plugin] BuildFieldHierarchy timeState=" << timeState
         << " meshCount=" << i.meshes.size() << "\n";
  struct LevelEntry {
    int level;
    std::string meshName;
  };

  std::unordered_map<std::string, std::vector<LevelEntry>> groupedMeshes;

  for (auto const &mesh_tuple : i.meshes) {
    const std::string openpmd_meshname = mesh_tuple.first;
    auto [baseName, level] = ParseMeshLevel(openpmd_meshname);
    groupedMeshes[baseName].push_back(LevelEntry{level, openpmd_meshname});
  }

  auto &hierarchyMap = meshHierarchyCache_.at(timeState);
  hierarchyMap.clear();

  for (auto &group : groupedMeshes) {
    debug2 << "[openpmd-api-plugin] Processing mesh group '" << group.first
           << "' levels=" << group.second.size() << "\n";
    std::vector<std::pair<int, std::string>> levels;
    levels.reserve(group.second.size());
    for (auto const &entry : group.second) {
      levels.emplace_back(entry.level, entry.meshName);
    }
    std::sort(levels.begin(), levels.end(),
              [](auto const &lhs, auto const &rhs) {
                return lhs.first < rhs.first;
              });

    std::string visitMeshName = group.first + "_mesh";
    MeshPatchHierarchy hierarchy =
        CreateHierarchyForGroup(visitMeshName, levels, i);

    hierarchy.metadataInitialized = true;
    hierarchyMap[visitMeshName] = hierarchy;

    meshMap_[visitMeshName] = std::tuple(DatasetType::Field, group.first);
    debug2 << "[openpmd-api-plugin] Registered AMR mesh '" << visitMeshName
           << "' backing openPMD mesh base '" << group.first
           << "' patches=" << hierarchy.patches.size()
           << " levels=" << hierarchy.numLevels << "\n";

    avtMeshMetaData *meshMd = new avtMeshMetaData;
    meshMd->name = visitMeshName;
    meshMd->meshType = AVT_AMR_MESH;
    meshMd->topologicalDimension = hierarchy.topologicalDim;
    meshMd->spatialDimension = hierarchy.spatialDim;
    meshMd->numBlocks = static_cast<int>(hierarchy.patches.size());
    meshMd->blockTitle = "patches";
    meshMd->blockPieceName = "patch";
    meshMd->numGroups = hierarchy.numLevels;
    meshMd->groupTitle = "levels";
    meshMd->groupPieceName = "level";
    meshMd->blockNames = hierarchy.blockNames;
    md->Add(meshMd);
    md->AddGroupInformation(hierarchy.numLevels,
                            static_cast<int>(hierarchy.patches.size()),
                            hierarchy.groupIds);

    // Register scalar variables using the finest-level mesh as representative.
    const std::string &representativeMeshName = levels.back().second;
    openPMD::Mesh mesh = i.meshes.at(representativeMeshName);

    if (mesh.scalar()) {
      std::string varname = group.first;
      avtCentering cent = GetCenteringType<openPMD::Mesh>(mesh);
      if (cent != AVT_UNKNOWN_CENT) {
        varMap_[varname] = std::tuple(visitMeshName,
                                       openPMD::MeshRecordComponent::SCALAR);
        AddScalarVarToMetaData(md, varname, visitMeshName, cent);
        debug2 << "[openpmd-api-plugin] Registered scalar var '" << varname
               << "' centering=" << CenteringToString(cent)
               << " representativeMesh='" << representativeMeshName << "'\n";
      } else {
        debug1 << "[openpmd-api-plugin] Skipping mesh '" << group.first
               << "' due to unsupported centering\n";
      }
    } else {
      for (auto const &rc : mesh) {
        const std::string &componentName = rc.first;
        std::string varname = group.first + "/" + componentName;
        avtCentering cent =
            GetCenteringType<openPMD::MeshRecordComponent>(rc.second);
        if (cent != AVT_UNKNOWN_CENT) {
          varMap_[varname] = std::tuple(visitMeshName, componentName);
          AddScalarVarToMetaData(md, varname, visitMeshName, cent);
          debug2 << "[openpmd-api-plugin] Registered component var '"
                 << varname << "' centering=" << CenteringToString(cent)
                 << " representativeMesh='" << representativeMeshName
                 << "'\n";
        } else {
          debug1 << "[openpmd-api-plugin] Skipping component '" << varname
                 << "' due to unsupported centering\n";
        }
      }
    }
  }

  debug1 << "[openpmd-api-plugin] Field hierarchy complete. Registered "
         << hierarchyMap.size() << " meshes and " << varMap_.size()
         << " scalars\n";
}

void avtopenpmdFileFormat::BuildParticleMetaData(avtDatabaseMetaData *md,
                                                 openPMD::Iteration const &i) {
  debug2 << "[openpmd-api-plugin] Skipping particle metadata population (not yet supported).\n";
  (void)md;
  (void)i;
}

// ****************************************************************************
//  Method: avtopenpmdFileFormat::PopulateDatabaseMetaData
//
//  Purpose:
//      This database meta-data object is like a table of contents for the
//      file.  By populating it, you are telling the rest of VisIt what
//      information it can request from you.
//
//  Programmer: benwibking -- generated by xml2avt
//  Creation:   Fri Dec 6 17:16:49 PST 2024
//
// ****************************************************************************

void avtopenpmdFileFormat::PopulateDatabaseMetaData(avtDatabaseMetaData *md,
                                                    int timeState) {
  debug1 << "[openpmd-api-plugin] PopulateDatabaseMetaData timeState="
         << timeState << "\n";
  // NOTE: openPMD's iteration index 'iter' can be an arbitrary number, and
  // can skip numbers. For instance, a dataset can have iterations = {550,
  // 600}.
  unsigned long long iterIdx = iterationIndex_.at(timeState);
  debug5 << "[openpmd-api-plugin] "
         << "PopulateDatabaseMetaData() for iteration " << iterIdx << "\n";

  // open openPMD::Iteration 'iter'
  openPMD::Iteration iter = series_.snapshots()[iterIdx];

  // populate field hierarchy metadata
  BuildFieldHierarchy(md, iter, timeState);

  // particles are currently unsupported in AMR mode
  BuildParticleMetaData(md, iter);

  debug1 << "[openpmd-api-plugin] PopulateDatabaseMetaData complete for"
         << " iteration=" << iterIdx << " meshes=" << meshMap_.size()
         << " vars=" << varMap_.size() << "\n";
}

void *avtopenpmdFileFormat::GetAuxiliaryData(const char *var, int timestep,
                                             int domain, const char *type,
                                             void *args, DestructorFunction &df) {
  (void)args;

  const char *varName = var != nullptr ? var : "<null>";
  const char *typeName = type != nullptr ? type : "<null>";
  debug1 << "[openpmd-api-plugin] GetAuxiliaryData type=" << typeName
         << " var=" << varName << " timestep=" << timestep
         << " domain=" << domain << "\n";

  if (type != nullptr &&
      (strcmp(type, AUXILIARY_DATA_DOMAIN_NESTING_INFORMATION) == 0 ||
       strcmp(type, AUXILIARY_DATA_DOMAIN_BOUNDARY_INFORMATION) == 0)) {
    if (var == nullptr) {
      debug1 << "[openpmd-api-plugin] Auxiliary request missing var name\n";
      return NULL;
    }

    if (timestep < 0 ||
        timestep >= static_cast<int>(meshHierarchyCache_.size())) {
      debug1 << "[openpmd-api-plugin] Auxiliary request out-of-range timestep\n";
      return NULL;
    }

    const std::string meshName(var);
    auto &hierarchyMap = meshHierarchyCache_[timestep];
    auto it = hierarchyMap.find(meshName);
    if (it == hierarchyMap.end()) {
      debug1 << "[openpmd-api-plugin] Auxiliary request missing mesh '"
             << meshName << "'\n";
      return NULL;
    }

    const MeshPatchHierarchy &hierarchy = it->second;
    if (hierarchy.patches.empty()) {
      debug1 << "[openpmd-api-plugin] Auxiliary request mesh '" << meshName
             << "' has no patches\n";
      return NULL;
    }

    if (strcmp(type, AUXILIARY_DATA_DOMAIN_NESTING_INFORMATION) == 0) {
      debug2 << "[openpmd-api-plugin] Building domain nesting for mesh '"
             << meshName << "'\n";
      avtStructuredDomainNesting *nesting = BuildDomainNesting(hierarchy);
      if (nesting == nullptr) {
        debug1 << "[openpmd-api-plugin] Domain nesting build returned nullptr\n";
        return NULL;
      }
      df = DeleteStructuredDomainNesting;
      debug1 << "[openpmd-api-plugin] Domain nesting ready for mesh '"
             << meshName << "'\n";
      return nesting;
    }

    if (strcmp(type, AUXILIARY_DATA_DOMAIN_BOUNDARY_INFORMATION) == 0) {
      debug2 << "[openpmd-api-plugin] Building domain boundaries for mesh '"
             << meshName << "'\n";
      avtStructuredDomainBoundaries *boundaries =
          BuildDomainBoundaries(hierarchy);
      if (boundaries == nullptr) {
        debug1 << "[openpmd-api-plugin] Domain boundary build returned nullptr\n";
        return NULL;
      }
      df = DeleteStructuredDomainBoundaries;
      debug1 << "[openpmd-api-plugin] Domain boundaries ready for mesh '"
             << meshName << "'\n";
      return boundaries;
    }
  }

  debug2 << "[openpmd-api-plugin] Forwarding auxiliary request to base class\n";
  return avtMTMDFileFormat::GetAuxiliaryData(var, timestep, domain, type, args,
                                             df);
}

std::pair<std::string, int>
avtopenpmdFileFormat::ParseMeshLevel(std::string const &meshName) const {
  int level = 0;
  std::string base = meshName;

  const std::string suffix = "_lvl";
  std::size_t pos = meshName.rfind(suffix);
  if (pos != std::string::npos) {
    std::size_t digitsBegin = pos + suffix.size();
    if (digitsBegin < meshName.size()) {
      std::size_t digitsEnd = digitsBegin;
      while (digitsEnd < meshName.size() && std::isdigit(meshName[digitsEnd])) {
        ++digitsEnd;
      }
      if (digitsEnd > digitsBegin) {
        level = std::stoi(meshName.substr(digitsBegin, digitsEnd - digitsBegin));
        base = meshName.substr(0, pos);
      }
    }
  }

  return {base, level};
}

MeshPatchHierarchy avtopenpmdFileFormat::CreateHierarchyForGroup(
    const std::string &visitMeshName,
    const std::vector<std::pair<int, std::string>> &levels,
    openPMD::Iteration const &iter) {
  MeshPatchHierarchy hierarchy;

  if (levels.empty()) {
    return hierarchy;
  }

  std::vector<int> uniqueLevels;
  uniqueLevels.reserve(levels.size());
  int lastLevel = std::numeric_limits<int>::min();
  for (auto const &levelPair : levels) {
    if (uniqueLevels.empty() || levelPair.first != lastLevel) {
      uniqueLevels.push_back(levelPair.first);
      lastLevel = levelPair.first;
    }
  }

  std::map<int, int> levelToGroupId;
  for (size_t idx = 0; idx < uniqueLevels.size(); ++idx) {
    levelToGroupId[uniqueLevels[idx]] = static_cast<int>(idx);
  }

  hierarchy.numLevels = static_cast<int>(uniqueLevels.size());
  hierarchy.levelValues = uniqueLevels;
  hierarchy.patchesPerLevel.assign(hierarchy.numLevels, {});
  hierarchy.levelCellSizes.assign(hierarchy.numLevels,
                                  std::array<double, 3>{0.0, 0.0, 0.0});

  std::map<int, int> patchCountByLevel;

  // Determine representative mesh to extract dimensional information.
  const std::string &representativeMeshName = levels.front().second;
  openPMD::Mesh representativeMesh = iter.meshes.at(representativeMeshName);
  GeometryData repGeom = GetGeometryXYZ(representativeMesh, true);
  int topoDim = 0;
  for (auto extent : repGeom.extent) {
    if (extent > 1) {
      ++topoDim;
    }
  }
  if (topoDim == 0) {
    topoDim = 1;
  }
  hierarchy.topologicalDim = topoDim;
  hierarchy.spatialDim = topoDim;

  debug5 << "[openpmd-api-plugin] BuildHierarchy -- visitMeshName="
         << visitMeshName << " levels=" << uniqueLevels.size()
         << " representativeMesh=" << representativeMeshName << "\n";
  debug5 << "[openpmd-api-plugin] Representative extent:";
  for (auto v : repGeom.extent) {
    debug5 << ' ' << v;
  }
  debug5 << " gridSpacing:";
  for (auto v : repGeom.gridSpacing) {
    debug5 << ' ' << v;
  }
  debug5 << " gridOrigin:";
  for (auto v : repGeom.gridOrigin) {
    debug5 << ' ' << v;
  }
  debug5 << "\n";

  auto getValueOr = [](const auto &container, int axis,
                       auto fallback) -> decltype(fallback) {
    if (axis >= 0 && axis < static_cast<int>(container.size())) {
      return container[axis];
    }
    return fallback;
  };

  for (auto const &[levelValue, meshName] : levels) {
    debug5 << "[openpmd-api-plugin]   processing mesh='" << meshName
           << "' level=" << levelValue << "\n";
    openPMD::Mesh mesh = iter.meshes.at(meshName);
    GeometryData geom = GetGeometryXYZ(mesh, true);
    double unitSI = mesh.gridUnitSI();
    avtCentering cent = GetCenteringType<openPMD::Mesh>(mesh);

    // Determine a representative record component for chunk enumeration.
    auto representativeComponent = mesh.scalar()
                                      ? mesh[openPMD::MeshRecordComponent::SCALAR]
                                      : mesh.begin()->second;

    openPMD::ChunkTable chunks = representativeComponent.availableChunks();
    debug5 << "[openpmd-api-plugin]     chunks available=" << chunks.size()
           << " unitSI=" << unitSI << " centering=" << cent << "\n";
    if (chunks.empty()) {
      openPMD::WrittenChunkInfo wholeDataset(
          openPMD::Offset(representativeComponent.getExtent().size(), 0),
          representativeComponent.getExtent());
      chunks.push_back(wholeDataset);
    }

    std::array<double, 3> spacingPerAxis{0.0, 0.0, 0.0};
    for (int axis = 0; axis < 3; ++axis) {
      double spacingValue = getValueOr(geom.gridSpacing, axis, 0.0);
      debug5 << "[openpmd-api-plugin]     axis=" << axis
             << " spacingRaw=" << spacingValue
             << " extent="
             << (axis < static_cast<int>(geom.extent.size())
                     ? geom.extent[axis]
                     : 0)
             << '\n';
      if (spacingValue == 0.0 &&
          axis < static_cast<int>(geom.extent.size()) &&
          geom.extent[axis] > 1) {
        debug5 << "[openpmd-api-plugin]       spacing fallback -> 1.0\n";
        spacingValue = 1.0;
      }
      spacingPerAxis[axis] = unitSI * spacingValue;
    }
    debug5 << "[openpmd-api-plugin]     spacingPerAxis=" << spacingPerAxis[0]
           << ',' << spacingPerAxis[1] << ',' << spacingPerAxis[2] << "\n";
    const int groupId = levelToGroupId[levelValue];
    hierarchy.levelCellSizes[groupId] = spacingPerAxis;

    for (auto const &chunk : chunks) {
      PatchInfo patch;
      patch.level = levelValue;
      patch.meshName = meshName;
      patch.offset = chunk.offset;
      patch.extent = chunk.extent;
      patch.centering = cent;

      for (int axis = 0; axis < 3; ++axis) {
        double spacing = spacingPerAxis[axis];
        patch.spacing[axis] = spacing;

        uint64_t offsetValue = 0;
        if (axis < static_cast<int>(patch.offset.size())) {
          offsetValue = patch.offset[axis];
        }

        double originValue = getValueOr(geom.gridOrigin, axis, 0.0);
        double gridOrigin = unitSI * originValue;
        patch.origin[axis] =
            gridOrigin + static_cast<double>(offsetValue) * spacing;
      }

      ComputeLogicalExtents(patch);

      hierarchy.patches.push_back(patch);
      const int patchIndex = static_cast<int>(hierarchy.patches.size()) - 1;
      std::ostringstream patchLabel;
      patchLabel << "Registered patch idx=" << patchIndex
                 << " levelValue=" << levelValue;
      LogPatchSummary(hierarchy.patches.back(), patchLabel.str());

      hierarchy.levelIdsPerPatch.push_back(groupId);
      hierarchy.groupIds.push_back(groupId);
      hierarchy.patchesPerLevel[groupId].push_back(patchIndex);

      int &patchCounter = patchCountByLevel[levelValue];
      std::ostringstream blockName;
      blockName << "level" << levelValue << ",patch" << patchCounter++;
      hierarchy.blockNames.push_back(blockName.str());
    }
  }

  hierarchy.levelRefinementRatios.clear();
  if (hierarchy.numLevels > 1) {
    for (int idx = 1; idx < hierarchy.numLevels; ++idx) {
      std::array<int, 3> ratio{1, 1, 1};
      const auto &coarseSpacing = hierarchy.levelCellSizes[idx - 1];
      const auto &fineSpacing = hierarchy.levelCellSizes[idx];
      for (int axis = 0; axis < 3; ++axis) {
        if (fineSpacing[axis] > 0.0 && coarseSpacing[axis] > 0.0) {
          int r = static_cast<int>(std::round(coarseSpacing[axis] /
                                              fineSpacing[axis]));
          ratio[axis] = std::max(1, r);
        }
      }
      hierarchy.levelRefinementRatios.push_back(ratio);
    }
  }

  debug5 << "[openpmd-api-plugin] Registered AMR mesh " << visitMeshName
         << " with " << hierarchy.patches.size() << " patches across "
         << hierarchy.numLevels << " levels.\n";

  debug2 << "[openpmd-api-plugin] Hierarchy summary mesh='" << visitMeshName
         << "' patches=" << hierarchy.patches.size()
         << " levels=" << hierarchy.numLevels << "\n";

  return hierarchy;
}

vtkDataSet *
avtopenpmdFileFormat::CreateRectilinearPatch(const PatchInfo &patch) const {
  vtkFloatArray *coords[3];
  int dimensions[3];

  for (int axis = 0; axis < 3; ++axis) {
    const uint64_t cells =
        axis < static_cast<int>(patch.extent.size()) ? patch.extent[axis] : 1;
    int nodes = static_cast<int>(cells);
    if (patch.centering == AVT_ZONECENT) {
      nodes = static_cast<int>(cells + 1);
    }
    if (nodes <= 0) {
      nodes = 1;
    }
    dimensions[axis] = nodes;

    coords[axis] = vtkFloatArray::New();
    coords[axis]->SetNumberOfTuples(nodes);
    float *array = static_cast<float *>(coords[axis]->GetVoidPointer(0));
    for (int idx = 0; idx < nodes; ++idx) {
      array[idx] = static_cast<float>(patch.origin[axis] +
                                      static_cast<double>(idx) *
                                          patch.spacing[axis]);
    }
  }

  vtkRectilinearGrid *grid = vtkRectilinearGrid::New();
  grid->SetDimensions(dimensions);
  grid->SetXCoordinates(coords[0]);
  grid->SetYCoordinates(coords[1]);
  grid->SetZCoordinates(coords[2]);

  for (int axis = 0; axis < 3; ++axis) {
    coords[axis]->Delete();
  }

  return grid;
}

avtStructuredDomainNesting *
avtopenpmdFileFormat::BuildDomainNesting(
    const MeshPatchHierarchy &hierarchy) const {
  if (hierarchy.patches.empty() || hierarchy.numLevels == 0) {
    debug1 << "[openpmd-api-plugin] BuildDomainNesting received empty hierarchy\n";
    return nullptr;
  }

  debug1 << "[openpmd-api-plugin] BuildDomainNesting patches="
         << hierarchy.patches.size() << " levels=" << hierarchy.numLevels
         << "\n";

  auto *nesting = new avtStructuredDomainNesting(
      static_cast<int>(hierarchy.patches.size()), hierarchy.numLevels);

  const int dims = std::max(1, hierarchy.spatialDim);
  nesting->SetNumDimensions(dims);

  nesting->SetLevelRefinementRatios(0, MakeRefinementVector({1, 1, 1}, dims));
  for (int level = 1; level < hierarchy.numLevels; ++level) {
    std::array<int, 3> ratio{1, 1, 1};
    if (level - 1 < static_cast<int>(hierarchy.levelRefinementRatios.size())) {
      ratio = hierarchy.levelRefinementRatios[level - 1];
    }
    nesting->SetLevelRefinementRatios(level, MakeRefinementVector(ratio, dims));
  }

  for (int level = 0; level < hierarchy.numLevels; ++level) {
    std::array<double, 3> sizes{0.0, 0.0, 0.0};
    if (level < static_cast<int>(hierarchy.levelCellSizes.size())) {
      sizes = hierarchy.levelCellSizes[level];
    }
    nesting->SetLevelCellSizes(level, MakeCellSizeVector(sizes, dims));
  }

  for (size_t patchIdx = 0; patchIdx < hierarchy.patches.size(); ++patchIdx) {
    const PatchInfo &patch = hierarchy.patches[patchIdx];
    const int levelIndex = hierarchy.levelIdsPerPatch[patchIdx];

    if (patchIdx == 0) {
      LogPatchSummary(patch, "BuildDomainNesting reference patch");
    }

    std::vector<int> children;
    if (levelIndex + 1 < hierarchy.numLevels) {
      const auto &candidateChildren =
          hierarchy.patchesPerLevel[levelIndex + 1];
      for (int childIdx : candidateChildren) {
        if (IsChildPatch(patch, hierarchy.patches[childIdx])) {
          children.push_back(childIdx);
        }
      }
    }

    std::vector<int> logicalExtents(6, 0);
    for (int axis = 0; axis < 3; ++axis) {
      logicalExtents[axis] = patch.logicalLower[axis];
      logicalExtents[axis + 3] = patch.logicalUpper[axis];
    }

    nesting->SetNestingForDomain(static_cast<int>(patchIdx), levelIndex,
                                 children, logicalExtents);
  }

  debug1 << "[openpmd-api-plugin] BuildDomainNesting completed\n";
  return nesting;
}

avtStructuredDomainBoundaries *
avtopenpmdFileFormat::BuildDomainBoundaries(
    const MeshPatchHierarchy &hierarchy) const {
  if (hierarchy.patches.empty()) {
    debug1 << "[openpmd-api-plugin] BuildDomainBoundaries received empty hierarchy\n";
    return nullptr;
  }

  debug1 << "[openpmd-api-plugin] BuildDomainBoundaries patches="
         << hierarchy.patches.size() << " levels=" << hierarchy.numLevels
         << "\n";

  bool canComputeNeighbors = true;
  auto *boundaries =
      new avtRectilinearDomainBoundaries(canComputeNeighbors);
  boundaries->SetNumDomains(static_cast<int>(hierarchy.patches.size()));

  if (hierarchy.numLevels > 1) {
    std::vector<int> refinement(hierarchy.numLevels - 1, 1);
    for (int level = 1; level < hierarchy.numLevels; ++level) {
      if (level - 1 < static_cast<int>(hierarchy.levelRefinementRatios.size())) {
        refinement[level - 1] = hierarchy.levelRefinementRatios[level - 1][0];
      }
    }
    boundaries->SetRefinementRatios(refinement);
  }

  for (size_t patchIdx = 0; patchIdx < hierarchy.patches.size(); ++patchIdx) {
    const PatchInfo &patch = hierarchy.patches[patchIdx];
    const int levelIndex = hierarchy.levelIdsPerPatch[patchIdx];

    if (patchIdx == 0) {
      LogPatchSummary(patch, "BuildDomainBoundaries reference patch");
    }

    int extents[6] = {0, 1, 0, 1, 0, 1};
    for (int axis = 0; axis < 3; ++axis) {
      if (axis < hierarchy.spatialDim) {
        extents[2 * axis] = patch.logicalLower[axis];
        extents[2 * axis + 1] = patch.logicalUpper[axis] + 1;
      }
    }

    boundaries->SetIndicesForAMRPatch(static_cast<int>(patchIdx), levelIndex,
                                      extents);
  }

  boundaries->CalculateBoundaries();
  debug1 << "[openpmd-api-plugin] BuildDomainBoundaries completed\n";
  return boundaries;
}

vtkDataArray *avtopenpmdFileFormat::LoadScalarPatchData(
    openPMD::Iteration const &iteration, const PatchInfo &patch,
    const std::string &component) const {
  debug1 << "[openpmd-api-plugin] LoadScalarPatchData mesh='" << patch.meshName
         << "' component='" << component << "' offset="
         << JoinContainer(patch.offset) << " extent="
         << JoinContainer(patch.extent) << "\n";

  auto meshIt = iteration.meshes.find(patch.meshName);
  if (meshIt == iteration.meshes.end()) {
    debug1 << "[openpmd-api-plugin] LoadScalarPatchData missing mesh '"
           << patch.meshName << "'\n";
    EXCEPTION1(InvalidVariableException, patch.meshName.c_str());
  }
  openPMD::Mesh mesh = meshIt->second;
  openPMD::MeshRecordComponent rcomp = mesh[component];

  uint64_t nelem = 1;
  for (auto const &nx : patch.extent) {
    nelem *= nx;
  }

  debug1 << "[openpmd-api-plugin] LoadScalarPatchData elements=" << nelem
         << " unitSI=" << rcomp.unitSI()
         << " datatype=" << static_cast<int>(rcomp.getDatatype()) << "\n";

  vtkDataArray *result = nullptr;

  if (rcomp.getDatatype() == openPMD::Datatype::DOUBLE) {
    vtkDoubleArray *data = vtkDoubleArray::New();
    data->SetNumberOfComponents(1);
    data->SetNumberOfTuples(nelem);
    double *buffer = data->GetPointer(0);
    try {
      rcomp.loadChunkRaw(buffer, patch.offset, patch.extent);
      const_cast<openPMD::Series &>(series_).flush();
    } catch (std::exception const &ex) {
      debug1 << "[openpmd-api-plugin] LoadScalarPatchData exception (double): "
             << ex.what() << "\n";
      data->Delete();
      throw;
    }
    ScaleVarData<double>(buffer, nelem, rcomp.unitSI());
    result = data;
  } else if (rcomp.getDatatype() == openPMD::Datatype::FLOAT) {
    vtkFloatArray *data = vtkFloatArray::New();
    data->SetNumberOfComponents(1);
    data->SetNumberOfTuples(nelem);
    float *buffer = data->GetPointer(0);
    try {
      rcomp.loadChunkRaw(buffer, patch.offset, patch.extent);
      const_cast<openPMD::Series &>(series_).flush();
    } catch (std::exception const &ex) {
      debug1 << "[openpmd-api-plugin] LoadScalarPatchData exception (float): "
             << ex.what() << "\n";
      data->Delete();
      throw;
    }
    ScaleVarData<float>(buffer, nelem, rcomp.unitSI());
    result = data;
  } else {
    debug1 << "[openpmd-api-plugin] LoadScalarPatchData unsupported datatype="
           << static_cast<int>(rcomp.getDatatype()) << "\n";
  }

  if (result != nullptr) {
    debug1 << "[openpmd-api-plugin] LoadScalarPatchData returning array\n";
  }

  return result;
}

vtkDataArray *avtopenpmdFileFormat::LoadVectorPatchData(
    openPMD::Iteration const &, const PatchInfo &,
    const std::vector<std::string> &) const {
  debug1 << "[openpmd-api-plugin] Vector variables are not supported\n";
  return nullptr;
}

vtkDataSet *avtopenpmdFileFormat::GetMesh(int timeState, int domain,
                                          const char *visit_meshname) {
  const char *meshName = visit_meshname != nullptr ? visit_meshname : "<null>";
  debug1 << "[openpmd-api-plugin] GetMesh timeState=" << timeState
         << " domain=" << domain << " mesh=" << meshName << "\n";

  if (visit_meshname == nullptr) {
    debug1 << "[openpmd-api-plugin] GetMesh received null mesh name\n";
    EXCEPTION1(InvalidVariableException, meshName);
  }

  auto hierarchyMapIt = meshHierarchyCache_.at(timeState).find(visit_meshname);
  if (hierarchyMapIt == meshHierarchyCache_.at(timeState).end()) {
    debug1 << "[openpmd-api-plugin] GetMesh missing mesh '" << meshName
           << "'\n";
    EXCEPTION1(InvalidVariableException, visit_meshname);
  }

  const MeshPatchHierarchy &hierarchy = hierarchyMapIt->second;
  if (domain < 0 || domain >= static_cast<int>(hierarchy.patches.size())) {
    debug1 << "[openpmd-api-plugin] GetMesh invalid domain index " << domain
           << " for mesh '" << meshName << "'\n";
    EXCEPTION1(InvalidVariableException, visit_meshname);
  }

  const PatchInfo &patch = hierarchy.patches.at(domain);
  debug5 << "[openpmd-api-plugin] GetMesh domain=" << domain
         << " level=" << patch.level << " using mesh " << patch.meshName
         << "\n";

  std::ostringstream ctx;
  ctx << "GetMesh returning patch domain=" << domain;
  LogPatchSummary(patch, ctx.str());

  vtkDataSet *grid = CreateRectilinearPatch(patch);
  debug1 << "[openpmd-api-plugin] GetMesh success mesh='" << meshName
         << "' domain=" << domain << "\n";
  return grid;
}

GeometryData avtopenpmdFileFormat::GetGeometry3D(openPMD::Mesh const &mesh,
                                                 bool insertMissingAxes) {
  // get dataOrder
  auto dataOrder = mesh.dataOrder();

  // get axis labels
  std::vector<std::string> axisLabels = mesh.axisLabels();
  debug2 << "[openpmd-api-plugin] GetGeometry3D axisLabelsRaw="
         << JoinStrings(axisLabels)
         << " dataOrder=" << DataOrderToString(dataOrder)
         << " insertMissingAxes=" << insertMissingAxes << "\n";
  if (doOverrideMeshAxisOrder_) {
    // for debugging
    debug5 << "Overriding Mesh axis labels. As-written axis labels: ";
    for (auto label : axisLabels) {
      debug5 << label << " ";
    }
    debug5 << '\n';

    axisLabels = overrideMeshAxisLabels_; // must be the same size as original
                                          // axisLabels!
    debug5 << "\tNew axis labels: ";
    for (auto label : axisLabels) {
      debug5 << label << " ";
    }
    debug5 << '\n';
  }

  // get array extents
  auto extent = mesh.getExtent();

  // get grid spacing
  std::vector<double> gridSpacing = mesh.gridSpacing<double>();
  // get grid origin
  std::vector<double> gridOrigin = mesh.gridGlobalOffset();

  // reverse ordering if mesh.DataOrder() == DataOrder::F
  if (dataOrder == openPMD::Mesh::DataOrder::F) {
    std::reverse(axisLabels.begin(), axisLabels.end());
    std::reverse(extent.begin(), extent.end());
    std::reverse(gridSpacing.begin(), gridSpacing.end());
    std::reverse(gridOrigin.begin(), gridOrigin.end());
  }

  if (insertMissingAxes) {
    // add any missing axis labels at the *beginning* with extent 1
    std::vector<std::string> canonicalAxes = {
        std::string("z"), std::string("y"), std::string("x")};
    for (auto const &axis : canonicalAxes) {
      if (std::find(axisLabels.begin(), axisLabels.end(), axis) ==
          axisLabels.end()) {
        axisLabels.insert(axisLabels.begin(), axis);
        extent.insert(extent.begin(), 1);
        gridSpacing.insert(gridSpacing.begin(), 0);
        gridOrigin.insert(gridOrigin.begin(), 0);
      }
    }
  }

  GeometryData geom;
  geom.axisLabels = axisLabels;
  geom.extent = extent;
  geom.gridOrigin = gridOrigin;
  geom.gridSpacing = gridSpacing;
  debug2 << "[openpmd-api-plugin] GetGeometry3D result axisLabels="
         << JoinStrings(geom.axisLabels) << " extent="
         << JoinContainer(geom.extent) << " spacing="
         << JoinContainer(geom.gridSpacing) << " origin="
         << JoinContainer(geom.gridOrigin) << "\n";
  return geom;
}

GeometryData avtopenpmdFileFormat::GetGeometryXYZ(openPMD::Mesh const &mesh,
                                                  bool insertMissingAxes) {
  GeometryData geom = GetGeometry3D(mesh, insertMissingAxes);

  // compute transposition of axisLabels -> {'x', 'y', 'z'}
  auto axisLabels = geom.axisLabels;
  debug5 << "[openpmd-api-plugin] Axis labels (slow->fast) before transpose: ";
  for (auto const &label : axisLabels) {
    debug5 << label << ' ';
  }
  debug5 << "\n";

  auto transpose = GetAxisTranspose(axisLabels, {"x", "y", "z"});
  debug5 << "[openpmd-api-plugin] Axis transpose mapping: ";
  for (auto idx : transpose) {
    debug5 << idx << ' ';
  }
  debug5 << "\n";

  // transpose geometry data
  TransposeVector(geom.axisLabels, transpose);
  TransposeVector(geom.extent, transpose);
  TransposeVector(geom.gridSpacing, transpose);
  TransposeVector(geom.gridOrigin, transpose);

  // verify we did it right
  assert(geom.axisLabels[0] == std::string("x"));
  assert(geom.axisLabels[1] == std::string("y"));
  assert(geom.axisLabels[2] == std::string("z"));

  debug2 << "[openpmd-api-plugin] GetGeometryXYZ result axisLabels="
         << JoinStrings(geom.axisLabels) << " extent="
         << JoinContainer(geom.extent) << " spacing="
         << JoinContainer(geom.gridSpacing) << " origin="
         << JoinContainer(geom.gridOrigin) << "\n";

  return geom;
}

// ****************************************************************************
//  Method: avtopenpmdFileFormat::GetVar
//
//  Purpose:
//      Gets a scalar variable associated with this file.  Although VTK has
//      support for many different types, the best bet is vtkFloatArray, since
//      that is supported everywhere through VisIt.
//
//  Arguments:
//      timestate  The index of the timestate.  If GetNTimesteps returned
//                 'N' time steps, this is guaranteed to be between 0 and N-1.
//      varname    The name of the variable requested.
//
//  Programmer: benwibking -- generated by xml2avt
//  Creation:   Fri Dec 6 17:16:49 PST 2024
//
// ****************************************************************************

vtkDataArray *avtopenpmdFileFormat::GetVar(int timeState, int domain,
                                           const char *varname) {
  const char *requestedVar = varname != nullptr ? varname : "<null>";
  debug1 << "[openpmd-api-plugin] GetVar timeState=" << timeState
         << " domain=" << domain << " var=" << requestedVar << "\n";
  if (varname == nullptr) {
    debug1 << "[openpmd-api-plugin] GetVar received null var name\n";
    EXCEPTION1(InvalidVariableException, requestedVar);
  }
  auto varIt = varMap_.find(varname);
  if (varIt == varMap_.end()) {
    debug1 << "[openpmd-api-plugin] GetVar missing var '" << requestedVar
           << "'\n";
    EXCEPTION1(InvalidVariableException, varname);
  }

  const std::string &visitMeshName = std::get<0>(varIt->second);
  const std::string &component = std::get<1>(varIt->second);

  auto &hierarchyMap = meshHierarchyCache_.at(timeState);
  auto meshIt = hierarchyMap.find(visitMeshName);
  if (meshIt == hierarchyMap.end()) {
    debug1 << "[openpmd-api-plugin] GetVar missing hierarchy for mesh '"
           << visitMeshName << "'\n";
    EXCEPTION1(InvalidVariableException, visitMeshName.c_str());
  }

  const MeshPatchHierarchy &hierarchy = meshIt->second;
  if (domain < 0 || domain >= static_cast<int>(hierarchy.patches.size())) {
    debug1 << "[openpmd-api-plugin] GetVar invalid domain index " << domain
           << " for mesh '" << visitMeshName << "'\n";
    EXCEPTION1(InvalidVariableException, visitMeshName.c_str());
  }

  const PatchInfo &patch = hierarchy.patches.at(domain);
  std::ostringstream ctx;
  ctx << "GetVar patch domain=" << domain << " component=" << component;
  LogPatchSummary(patch, ctx.str());

  unsigned long long iter = iterationIndex_.at(timeState);
  debug1 << "[openpmd-api-plugin] GetVar loading iteration=" << iter
         << " component='" << component << "'\n";
  openPMD::Iteration iteration = series_.snapshots()[iter];

  vtkDataArray *data = LoadScalarPatchData(iteration, patch, component);
  if (data == nullptr) {
    debug1 << "[openpmd-api-plugin] GetVar LoadScalarPatchData returned nullptr"
           << " for var '" << requestedVar << "'\n";
    EXCEPTION1(InvalidVariableException, varname);
  }

  debug1 << "[openpmd-api-plugin] GetVar success var='" << requestedVar
         << "'\n";

  return data;
}

// ****************************************************************************
//  Method: avtopenpmdFileFormat::GetVectorVar
//
//  Purpose:
//      Gets a vector variable associated with this file.  Although VTK has
//      support for many different types, the best bet is vtkFloatArray, since
//      that is supported everywhere through VisIt.
//
//  Arguments:
//      timestate  The index of the timestate.  If GetNTimesteps returned
//                 'N' time steps, this is guaranteed to be between 0 and N-1.
//      varname    The name of the variable requested.
//
//  Programmer: benwibking -- generated by xml2avt
//  Creation:   Fri Dec 6 17:16:49 PST 2024
//
// ****************************************************************************

vtkDataArray *avtopenpmdFileFormat::GetVectorVar(int, int,
                                                 const char *varname) {
  EXCEPTION1(InvalidVariableException, varname);
}

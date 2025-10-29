// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ****************************************************************************
//  avtopenpmdFileFormat.C
// ****************************************************************************

#include <algorithm>
#include <array>
#include <cctype>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <csignal>
#include <atomic>
#include <iostream>
#include <map>
#include <numeric>
#include <limits>
#include <sstream>
#include <cstring>
#include <unordered_map>
#include <stdexcept>
#include <string>
#include <unistd.h>

#include <avtDatabaseMetaData.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkCellType.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkRectilinearGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnsignedCharArray.h>

#include <openPMD/openPMD.hpp>

#include <cpptrace/cpptrace.hpp>

#include <DBOptionsAttributes.h>
#include <DebugStream.h>
#include <Expression.h>
#include <InvalidVariableException.h>
#include <InvalidFilesException.h>
#include <avtGhostData.h>

#include "avtopenpmdFileFormat.h"

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

  bool fineStrictlyInside = true;
  bool boxesOverlap = true;
  bool hasFinerResolutionAxis = false;

  for (int axis = 0; axis < 3; ++axis) {
    const double coarseSpacing = std::abs(coarse.spacing[axis]);
    const double fineSpacing = std::abs(fine.spacing[axis]);
    const double tol = std::max({coarseSpacing, fineSpacing, 1.0}) * 1e-5;

    if (fineMin[axis] < coarseMin[axis] - tol ||
        fineMax[axis] > coarseMax[axis] + tol) {
      fineStrictlyInside = false;
    }

    const double coarseWidth = coarseMax[axis] - coarseMin[axis];
    const double fineWidth = fineMax[axis] - fineMin[axis];
    bool axisOverlap = false;
    if (coarseWidth <= tol && fineWidth <= tol) {
      const double minDelta = std::abs(fineMin[axis] - coarseMin[axis]);
      const double maxDelta = std::abs(fineMax[axis] - coarseMax[axis]);
      axisOverlap = (minDelta <= tol) && (maxDelta <= tol);
    } else {
      axisOverlap = (fineMin[axis] < coarseMax[axis] - tol) &&
                    (fineMax[axis] > coarseMin[axis] + tol);
    }
    if (!axisOverlap) {
      boxesOverlap = false;
    }

    if (coarseSpacing > 0.0 && fineSpacing > 0.0 &&
        fineSpacing + tol < coarseSpacing) {
      hasFinerResolutionAxis = true;
    }
  }

  if (!boxesOverlap) {
    return false;
  }

  if (fineStrictlyInside) {
    return true;
  }

  if (!hasFinerResolutionAxis) {
    return false;
  }

  return true;
}

inline std::pair<std::string, std::string>
SplitRecordComponentPath(const std::string &path) {
  auto delimiter = path.find('/');
  if (delimiter == std::string::npos) {
    return {path, openPMD::RecordComponent::SCALAR};
  }

  std::string record = path.substr(0, delimiter);
  std::string component = path.substr(delimiter + 1);
  if (component.empty()) {
    component = openPMD::RecordComponent::SCALAR;
  }
  return {record, component};
}

inline std::size_t ElementCountFromExtent(const openPMD::Extent &extent) {
  if (extent.empty()) {
    return 0;
  }
  std::size_t count = 1;
  for (auto dim : extent) {
    count *= static_cast<std::size_t>(dim == 0 ? 1 : dim);
  }
  return count;
}

template <typename T>
inline void DuplicateHighEndNodes(const PatchInfo &patch, T *values) {
  if (values == nullptr || patch.centering != AVT_NODECENT) {
    return;
  }

  std::array<uint64_t, 3> vtkDims{1, 1, 1};
  for (int axis = 0; axis < 3; ++axis) {
    if (axis < static_cast<int>(patch.extent.size())) {
      uint64_t dim = patch.extent[axis];
      vtkDims[axis] = dim == 0 ? 1 : dim;
    }
  }

  const auto &storageDims = patch.storageExtentCanonical;
  const auto &storageToVtk = patch.storageToVtk;
  std::array<uint64_t, 3> validDims{1, 1, 1};
  bool requiresDuplication = false;
  for (int axis = 0; axis < 3; ++axis) {
    int storageIndex = storageToVtk[axis];
    uint64_t dimValue = vtkDims[axis];
    if (storageIndex >= 0 &&
        storageIndex < static_cast<int>(storageDims.size())) {
      dimValue = storageDims[static_cast<size_t>(storageIndex)];
    }
    if (dimValue == 0) {
      dimValue = 1;
    }
    validDims[axis] = dimValue;
    if (vtkDims[axis] > dimValue) {
      requiresDuplication = true;
    }
  }

  if (!requiresDuplication) {
    return;
  }

  uint64_t totalElements = vtkDims[0] * vtkDims[1] * vtkDims[2];
  for (uint64_t linear = 0; linear < totalElements; ++linear) {
    uint64_t tmp = linear;
    std::array<uint64_t, 3> coords{0, 0, 0};
    for (int axis = 0; axis < 3; ++axis) {
      const uint64_t dim = vtkDims[axis];
      if (dim > 0) {
        coords[axis] = tmp % dim;
        tmp /= dim;
      } else {
        coords[axis] = 0;
      }
    }

    bool needsCopy = false;
    std::array<uint64_t, 3> srcCoords = coords;
    for (int axis = 0; axis < 3; ++axis) {
      uint64_t valid = validDims[axis];
      if (coords[axis] >= valid) {
        needsCopy = true;
        srcCoords[axis] = valid > 0 ? valid - 1 : 0;
      }
    }

    if (!needsCopy) {
      continue;
    }

    uint64_t srcIndex =
        srcCoords[0] + vtkDims[0] * (srcCoords[1] + vtkDims[1] * srcCoords[2]);
    values[linear] = values[static_cast<size_t>(srcIndex)];
  }
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

inline bool IsAbsolutePath(const std::string &path) {
  if (path.empty()) {
    return false;
  }
#ifdef _WIN32
  if (path.size() >= 2 &&
      std::isalpha(static_cast<unsigned char>(path[0])) &&
      path[1] == ':') {
    return true;
  }
  if (path.size() >= 2 && path[0] == '\\' && path[1] == '\\') {
    return true;
  }
  return path[0] == '/' || path[0] == '\\';
#else
  return path[0] == '/';
#endif
}

inline std::string ParentDirectory(const std::string &path) {
  if (path.empty()) {
    return "";
  }

  const auto pos = path.find_last_of("/\\");
  if (pos == std::string::npos) {
    return "";
  }

#ifdef _WIN32
  if (pos == 2 && path[1] == ':') {
    return path.substr(0, pos + 1);
  }
#endif
  if (pos == 0) {
    return path.substr(0, 1);
  }
  return path.substr(0, pos);
}

inline bool EndsWithSeparator(const std::string &path) {
  if (path.empty()) {
    return false;
  }
  const char last = path.back();
  return last == '/' || last == '\\';
}

inline std::string JoinPath(const std::string &parent,
                            const std::string &child) {
  if (child.empty()) {
    return parent;
  }
  if (IsAbsolutePath(child) || parent.empty()) {
    return child;
  }

  std::string result = parent;
  if (!EndsWithSeparator(result)) {
#ifdef _WIN32
    if (result.size() == 2 && result[1] == ':') {
      result.push_back('\\');
    } else {
      result.push_back('/');
    }
#else
    result.push_back('/');
#endif
  }
  result += child;
  return result;
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
  std::string opmd_overrideMeshDataOrder;

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

    if (std::getline(file, opmd_overrideMeshDataOrder)) {
      std::istringstream orderStream(opmd_overrideMeshDataOrder);
      std::string token;
      if (orderStream >> token) {
        if (!token.empty() && token[0] != '#') {
          char letter = static_cast<char>(std::toupper(static_cast<unsigned char>(token[0])));
          if (letter == 'F') {
            doOverrideMeshDataOrder_ = true;
            overrideMeshDataOrder_ = openPMD::Mesh::DataOrder::F;
            debug3 << "[openpmd-api-plugin] overrideMeshDataOrder enabled: F\n";
          } else if (letter == 'C') {
            doOverrideMeshDataOrder_ = true;
            overrideMeshDataOrder_ = openPMD::Mesh::DataOrder::C;
            debug3 << "[openpmd-api-plugin] overrideMeshDataOrder enabled: C\n";
          } else {
            debug1 << "[openpmd-api-plugin] Unknown dataOrder override token '"
                   << token << "'\n";
          }
        } else {
          debug3 << "[openpmd-api-plugin] overrideMeshDataOrder ignored (comment)\n";
        }
      } else {
        debug3 << "[openpmd-api-plugin] overrideMeshDataOrder disabled\n";
      }
    } else {
      debug3 << "[openpmd-api-plugin] overrideMeshDataOrder disabled\n";
    }

    // close file
    file.close();
  }

  // construct complete filepath
  std::string descriptor_path(filename);
  std::string parent_path = ParentDirectory(descriptor_path);
  std::string opmd_filepath = JoinPath(parent_path, opmd_filestring);
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
  particleCache_.resize(numTimesteps);

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
  for (auto &timestepCache : particleCache_) {
    for (auto &speciesEntry : timestepCache) {
      speciesEntry.second.cachedMesh = nullptr;
      speciesEntry.second.meshCached = false;
    }
  }
}

void avtopenpmdFileFormat::PopulateHierarchyCache(
    openPMD::Iteration const &i, int timeState, avtDatabaseMetaData *md) {
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

    if (md != nullptr) {
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
      meshMd->containsGhostZones = AVT_HAS_GHOSTS;
      meshMd->presentGhostZoneTypes = (1 << AVT_NESTING_GHOST_ZONES);
      md->Add(meshMd);
      md->AddGroupInformation(hierarchy.numLevels,
                              static_cast<int>(hierarchy.patches.size()),
                              hierarchy.groupIds);
    }

    const std::string &representativeMeshName = levels.back().second;
    openPMD::Mesh mesh = i.meshes.at(representativeMeshName);

    if (mesh.scalar()) {
      std::string varname = group.first;
      avtCentering cent = GetCenteringType<openPMD::Mesh>(mesh);
      if (cent != AVT_UNKNOWN_CENT) {
        varMap_[varname] = std::tuple(visitMeshName,
                                       openPMD::MeshRecordComponent::SCALAR);
        if (md != nullptr) {
          AddScalarVarToMetaData(md, varname, visitMeshName, cent);
        }
        debug2 << "[openpmd-api-plugin] Registered scalar var '" << varname
               << "' centering=" << CenteringToString(cent)
               << " representativeMesh='" << representativeMeshName << "'\n";
      } else {
        debug1 << "[openpmd-api-plugin] Skipping mesh '" << group.first
               << "' due to unsupported centering\n";
      }
    } else {
      std::vector<std::pair<std::string, openPMD::MeshRecordComponent>>
          components;
      components.reserve(mesh.size());
      for (auto const &rc : mesh) {
        components.emplace_back(rc.first, rc.second);
      }

      if (components.empty()) {
        debug1 << "[openpmd-api-plugin] Mesh '" << group.first
               << "' has no record components; skipping registration\n";
      } else {
        const size_t componentCount = components.size();
        avtCentering vectorCentering = GetCenteringType<openPMD::Mesh>(mesh);
        if (vectorCentering == AVT_UNKNOWN_CENT) {
          vectorCentering =
              GetCenteringType<openPMD::MeshRecordComponent>(components.front().second);
        }

        if (vectorCentering == AVT_UNKNOWN_CENT) {
          debug1 << "[openpmd-api-plugin] Skipping vector var '" << group.first
                 << "' due to unsupported centering\n";
        } else if (componentCount ==
                   static_cast<size_t>(hierarchy.spatialDim)) {
          std::vector<std::string> componentNames;
          componentNames.reserve(componentCount);
          for (auto const &entry : components) {
            componentNames.push_back(entry.first);
          }

          vectorVarMap_[group.first] =
              std::make_tuple(visitMeshName, componentNames);
          if (md != nullptr) {
            AddVectorVarToMetaData(md, group.first, visitMeshName,
                                   vectorCentering,
                                   static_cast<int>(componentCount));
          }
          debug2 << "[openpmd-api-plugin] Registered vector var '"
                 << group.first << "' components="
                 << JoinStrings(componentNames)
                 << " centering=" << CenteringToString(vectorCentering)
                 << " representativeMesh='" << representativeMeshName
                 << "'\n";
        } else {
          debug1 << "[openpmd-api-plugin] Mesh '" << group.first
                 << "' component count " << componentCount
                 << " does not match spatial dimension "
                 << hierarchy.spatialDim << "; registering components only\n";
        }

        for (auto const &entry : components) {
          const std::string &componentName = entry.first;
          std::string varname = group.first + "/" + componentName;
          avtCentering cent =
              GetCenteringType<openPMD::MeshRecordComponent>(entry.second);
          if (cent != AVT_UNKNOWN_CENT) {
            varMap_[varname] = std::tuple(visitMeshName, componentName);
            if (md != nullptr) {
              AddScalarVarToMetaData(md, varname, visitMeshName, cent);
            }
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
  }

  if (md == nullptr) {
    debug1 << "[openpmd-api-plugin] Field hierarchy cache populated for"
           << " timeState=" << timeState
           << " meshes=" << hierarchyMap.size() << "\n";
  }
}

void avtopenpmdFileFormat::BuildFieldHierarchy(avtDatabaseMetaData *md,
                                               openPMD::Iteration const &i,
                                               int timeState) {
  debug1 << "[openpmd-api-plugin] BuildFieldHierarchy timeState=" << timeState
         << " meshCount=" << i.meshes.size() << "\n";

  PopulateHierarchyCache(i, timeState, md);

  auto &hierarchyMap = meshHierarchyCache_.at(timeState);
  debug1 << "[openpmd-api-plugin] Field hierarchy complete. Registered "
         << hierarchyMap.size() << " meshes, " << varMap_.size()
         << " scalar components and " << vectorVarMap_.size()
         << " vector variables\n";
}

void avtopenpmdFileFormat::BuildParticleMetaData(avtDatabaseMetaData *md,
                                                 openPMD::Iteration const &i,
                                                 int timeState) {
  debug1 << "[openpmd-api-plugin] BuildParticleMetaData timeState=" << timeState
         << " speciesCount=" << i.particles.size() << "\n";

  if (timeState < 0 ||
      timeState >= static_cast<int>(particleCache_.size())) {
    debug1 << "[openpmd-api-plugin] BuildParticleMetaData invalid timeState="
           << timeState << " cacheSize=" << particleCache_.size() << "\n";
    return;
  }

  auto &speciesCache = particleCache_.at(timeState);
  speciesCache.clear();

  for (auto const &speciesEntry : i.particles) {
    const std::string &speciesName = speciesEntry.first;
    const std::string visitMeshName = speciesName + "_particles";
    const openPMD::ParticleSpecies &species = speciesEntry.second;

    ParticleSpeciesInfo info;
    info.speciesName = speciesName;
    info.visitMeshName = visitMeshName;

    auto positionIt = species.find("position");
    if (positionIt != species.end()) {
      const openPMD::Record &positionRecord = positionIt->second;

      info.positionComponents.clear();
      info.positionComponents.reserve(positionRecord.size());
      for (auto const &componentEntry : positionRecord) {
        info.positionComponents.push_back(componentEntry.first);
      }
      if (info.positionComponents.empty() && positionRecord.scalar()) {
        info.positionComponents.push_back(openPMD::RecordComponent::SCALAR);
      }

      info.storageAxisLabels = info.positionComponents;
      info.storageToVtk.assign(info.positionComponents.size(), -1);

      std::vector<std::string> desiredAxisOrder;
      if (doOverrideParticleAxisOrder_ && !overrideParticleAxisLabels_.empty()) {
        desiredAxisOrder = overrideParticleAxisLabels_;
      } else {
        desiredAxisOrder = {std::string("x"), std::string("y"),
                            std::string("z")};
      }

      for (std::size_t storageIdx = 0;
           storageIdx < info.positionComponents.size(); ++storageIdx) {
        const std::string &componentName = info.positionComponents[storageIdx];

        int mappedAxis = -1;
        if (!desiredAxisOrder.empty()) {
          auto desiredIt = std::find(desiredAxisOrder.begin(),
                                     desiredAxisOrder.end(), componentName);
          if (desiredIt != desiredAxisOrder.end()) {
            mappedAxis =
                static_cast<int>(std::distance(desiredAxisOrder.begin(), desiredIt));
          }
        }

        if (!doOverrideParticleAxisOrder_ && mappedAxis == -1) {
          std::string lowered = componentName;
          std::transform(lowered.begin(), lowered.end(), lowered.begin(),
                         [](unsigned char c) { return std::tolower(c); });
          if (lowered == "x") {
            mappedAxis = 0;
          } else if (lowered == "y") {
            mappedAxis = 1;
          } else if (lowered == "z") {
            mappedAxis = 2;
          }
        }

        if (mappedAxis >= 0 && mappedAxis < 3) {
          info.storageToVtk[storageIdx] = mappedAxis;
        }
      }

      info.spatialDim = 0;
      for (int axis : info.storageToVtk) {
        if (axis >= 0 && axis < 3) {
          info.spatialDim = std::max(info.spatialDim, axis + 1);
        }
      }
      if (info.spatialDim == 0) {
        info.spatialDim = static_cast<int>(
            std::min<std::size_t>(3, info.positionComponents.size()));
      }

      if (!info.positionComponents.empty()) {
        std::size_t particleCount = 0;
        if (info.positionComponents.front() == openPMD::RecordComponent::SCALAR) {
          particleCount =
              ElementCountFromExtent(positionRecord.getExtent());
        } else {
        const std::string &componentName = info.positionComponents.front();
        auto compIt = positionRecord.find(componentName);
        if (compIt != positionRecord.end()) {
          particleCount = ElementCountFromExtent(compIt->second.getExtent());
        }
        }
        info.particleCount = particleCount;
      }
    } else {
      debug2 << "[openpmd-api-plugin] Particle species '" << speciesName
             << "' missing position record; treating as empty.\n";
      info.storageToVtk.clear();
      info.spatialDim = 0;
      info.particleCount = 0;
    }

    info.metadataInitialized = true;
    speciesCache[visitMeshName] = info;

    meshMap_[visitMeshName] =
        std::tuple(DatasetType::ParticleSpecies, speciesName);

    if (md != nullptr) {
      avtMeshMetaData *meshMd = new avtMeshMetaData;
      meshMd->name = visitMeshName;
      meshMd->meshType = AVT_POINT_MESH;
      meshMd->topologicalDimension = 0;
      meshMd->spatialDimension = std::max(1, info.spatialDim);
      meshMd->numBlocks = 1;
      meshMd->blockTitle = "species";
      meshMd->blockPieceName = "species";
      meshMd->numGroups = 0;
      meshMd->blockNames.push_back(speciesName);
      meshMd->containsGhostZones = AVT_NO_GHOSTS;
      meshMd->presentGhostZoneTypes = 0;
      md->Add(meshMd);
    }

    for (auto const &recordEntry : species) {
      const std::string &recordName = recordEntry.first;
      if (recordName == "position") {
        continue;
      }
      const openPMD::Record &record = recordEntry.second;

      if (record.scalar()) {
        std::string componentPath =
            recordName + "/" + openPMD::RecordComponent::SCALAR;
        std::string varName = speciesName + "/" + recordName;
        varMap_[varName] = std::tuple(visitMeshName, componentPath);
        if (md != nullptr) {
          AddScalarVarToMetaData(md, varName, visitMeshName, AVT_NODECENT);
        }
        debug2 << "[openpmd-api-plugin] Registered particle scalar '"
               << varName << "'\n";
      } else {
        std::vector<std::string> componentPaths;
        componentPaths.reserve(record.size());
        for (auto const &componentEntry : record) {
          const std::string &componentName = componentEntry.first;
          std::string componentPath = recordName + "/" + componentName;
          std::string scalarName = speciesName + "/" + recordName + "/" +
                                   componentName;
          componentPaths.push_back(componentPath);
          varMap_[scalarName] = std::tuple(visitMeshName, componentPath);
          if (md != nullptr) {
            AddScalarVarToMetaData(md, scalarName, visitMeshName,
                                   AVT_NODECENT);
          }
          debug2 << "[openpmd-api-plugin] Registered particle component '"
                 << scalarName << "'\n";
        }

        if (!componentPaths.empty()) {
          std::string vectorName = speciesName + "/" + recordName;
          vectorVarMap_[vectorName] =
              std::tuple(visitMeshName, componentPaths);
          if (md != nullptr) {
            AddVectorVarToMetaData(md, vectorName, visitMeshName,
                                   AVT_NODECENT,
                                   static_cast<int>(componentPaths.size()));
          }
          debug2 << "[openpmd-api-plugin] Registered particle vector '"
                 << vectorName << "' components="
                 << JoinStrings(componentPaths) << "\n";
        }
      }
    }
  }

  debug1 << "[openpmd-api-plugin] Particle metadata complete timeState="
         << timeState << " meshes=" << speciesCache.size() << "\n";
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
  BuildParticleMetaData(md, iter, timeState);

  debug1 << "[openpmd-api-plugin] PopulateDatabaseMetaData complete for"
         << " iteration=" << iterIdx << " meshes=" << meshMap_.size()
         << " vars=" << varMap_.size() << "\n";
}

void avtopenpmdFileFormat::EnsureHierarchyInitialized(int timeState) {
  if (timeState < 0 ||
      timeState >= static_cast<int>(meshHierarchyCache_.size())) {
    debug1 << "[openpmd-api-plugin] EnsureHierarchyInitialized out-of-range"
           << " timeState=" << timeState
           << " cacheSize=" << meshHierarchyCache_.size() << "\n";
    EXCEPTION1(InvalidVariableException, "timestep out of range");
  }

  auto &hierarchyMap = meshHierarchyCache_.at(timeState);
  if (!hierarchyMap.empty()) {
    return;
  }

  unsigned long long iterIdx = iterationIndex_.at(timeState);
  debug1 << "[openpmd-api-plugin] Lazy hierarchy initialization for timeState="
         << timeState << " iteration=" << iterIdx << "\n";
  openPMD::Iteration iter = series_.snapshots()[iterIdx];
  PopulateHierarchyCache(iter, timeState, nullptr);
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

    EnsureHierarchyInitialized(timestep);

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
      if (domain == -1) {
        debug1 << "[openpmd-api-plugin] Global boundary request not supported, returning null\n";
        return NULL;
      }
      if (domain < 0 ||
          domain >= static_cast<int>(hierarchy.patches.size())) {
        debug1 << "[openpmd-api-plugin] Auxiliary request invalid domain "
               << domain << " for mesh '" << meshName << "'\n";
        return NULL;
      }
      debug2 << "[openpmd-api-plugin] Building domain boundary list for mesh '"
             << meshName << "' domain=" << domain << "\n";
      avtLocalStructuredDomainBoundaryList *boundaryList =
          BuildDomainBoundaryList(hierarchy, domain);
      if (boundaryList == nullptr) {
        debug1 << "[openpmd-api-plugin] Domain boundary build returned nullptr\n";
        return NULL;
      }
      df = avtLocalStructuredDomainBoundaryList::Destruct;
      debug1 << "[openpmd-api-plugin] Domain boundary list ready for mesh '"
             << meshName << "' domain=" << domain << "\n";
      return boundaryList;
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

  auto getValueOr = [](const auto &container, int axis,
                       auto fallback) -> decltype(fallback) {
    if (axis >= 0 && axis < static_cast<int>(container.size())) {
      return container[axis];
    }
    return fallback;
  };

  std::vector<int> uniqueLevels;
  uniqueLevels.reserve(levels.size());
  std::map<int, std::array<double, 3>> levelSpacingMap;
  std::map<int, std::string> levelToMeshName;

  for (auto const &levelPair : levels) {
    const int levelValue = levelPair.first;
    const std::string &meshName = levelPair.second;
    if (levelSpacingMap.find(levelValue) != levelSpacingMap.end()) {
      continue;
    }

    openPMD::Mesh spacingMesh = iter.meshes.at(meshName);
    openPMD::MeshRecordComponent representativeComponentForSpacing =
        spacingMesh.scalar()
            ? spacingMesh[openPMD::MeshRecordComponent::SCALAR]
            : spacingMesh.begin()->second;
    GeometryData spacingGeom =
        GetGeometryXYZ(spacingMesh, &representativeComponentForSpacing, true);

    const double unitSI = spacingMesh.gridUnitSI();
    std::array<double, 3> spacingPerAxis{0.0, 0.0, 0.0};
    for (int axis = 0; axis < 3; ++axis) {
      double spacingValue = getValueOr(spacingGeom.gridSpacing, axis, 0.0);
      if (spacingValue == 0.0 &&
          axis < static_cast<int>(spacingGeom.extent.size()) &&
          spacingGeom.extent[axis] > 1) {
        spacingValue = 1.0;
      }
      spacingPerAxis[axis] = unitSI * spacingValue;
    }

    uniqueLevels.push_back(levelValue);
    levelSpacingMap[levelValue] = spacingPerAxis;
    levelToMeshName[levelValue] = meshName;
  }

  if (uniqueLevels.empty()) {
    return hierarchy;
  }

  auto spacingMetric = [](const std::array<double, 3> &spacing) {
    double metric = 0.0;
    for (double value : spacing) {
      double absValue = std::abs(value);
      if (absValue > metric) {
        metric = absValue;
      }
    }
    return metric;
  };

  std::vector<int> sortedLevels = uniqueLevels;
  std::sort(sortedLevels.begin(), sortedLevels.end(), [&](int lhs, int rhs) {
    const double lhsMetric = spacingMetric(levelSpacingMap[lhs]);
    const double rhsMetric = spacingMetric(levelSpacingMap[rhs]);
    if (std::abs(lhsMetric - rhsMetric) > 1e-12) {
      return lhsMetric > rhsMetric;
    }
    return lhs < rhs;
  });

  std::map<int, int> levelToGroupId;
  for (size_t idx = 0; idx < sortedLevels.size(); ++idx) {
    levelToGroupId[sortedLevels[idx]] = static_cast<int>(idx);
  }

  hierarchy.numLevels = static_cast<int>(sortedLevels.size());
  hierarchy.levelValues = sortedLevels;
  hierarchy.patchesPerLevel.assign(hierarchy.numLevels, {});
  hierarchy.levelCellSizes.assign(hierarchy.numLevels,
                                  std::array<double, 3>{0.0, 0.0, 0.0});

  std::map<int, int> patchCountByLevel;

  // Determine representative mesh to extract dimensional information.
  const std::string &representativeMeshName =
      levelToMeshName.at(sortedLevels.front());
  openPMD::Mesh representativeMesh = iter.meshes.at(representativeMeshName);
  openPMD::MeshRecordComponent representativeComponentForGeom =
      representativeMesh.scalar()
          ? representativeMesh[openPMD::MeshRecordComponent::SCALAR]
          : representativeMesh.begin()->second;
  GeometryData repGeom =
      GetGeometryXYZ(representativeMesh, &representativeComponentForGeom, true);
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
         << visitMeshName << " levels=" << hierarchy.numLevels
         << " levelOrder=" << JoinContainer(hierarchy.levelValues)
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

  for (auto const &[levelValue, meshName] : levels) {
    debug5 << "[openpmd-api-plugin]   processing mesh='" << meshName
           << "' level=" << levelValue << "\n";
    openPMD::Mesh mesh = iter.meshes.at(meshName);
    // Determine a representative record component for chunk enumeration.
    auto representativeComponent = mesh.scalar()
                                      ? mesh[openPMD::MeshRecordComponent::SCALAR]
                                      : mesh.begin()->second;

    GeometryData geom = GetGeometryXYZ(mesh, &representativeComponent, true);
    auto dataOrder = GetEffectiveMeshDataOrder(mesh);
    double unitSI = mesh.gridUnitSI();
    avtCentering cent = GetCenteringType<openPMD::Mesh>(mesh);

    std::array<uint64_t, 3> globalExtent{{1, 1, 1}};
    for (int axis = 0; axis < 3; ++axis) {
      if (axis < static_cast<int>(geom.extent.size())) {
        uint64_t extentValue = geom.extent[axis];
        if (extentValue == 0) {
          extentValue = 1;
        }
        globalExtent[axis] = extentValue;
      }
    }

    std::vector<std::string> rawAxisLabels = mesh.axisLabels();
    if (dataOrder == openPMD::Mesh::DataOrder::F) {
      std::reverse(rawAxisLabels.begin(), rawAxisLabels.end());
    }

    const std::vector<std::string> &storageAxisLabels =
        geom.storageAxisLabels.empty() ? rawAxisLabels : geom.storageAxisLabels;
    std::vector<int> storageToVtkVec = geom.storageToVtk;
    if (storageToVtkVec.empty()) {
      storageToVtkVec.resize(3);
      std::iota(storageToVtkVec.begin(), storageToVtkVec.end(), 0);
    }

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

    avtCentering chunkCentering = cent;
    if (chunkCentering == AVT_UNKNOWN_CENT) {
      chunkCentering =
          GetCenteringType<openPMD::MeshRecordComponent>(representativeComponent);
    }

    for (auto const &chunk : chunks) {
      PatchInfo patch;
      patch.level = levelValue;
      patch.meshName = meshName;
      patch.centering = chunkCentering;
      patch.dataOrder = dataOrder;
      patch.storageAxisLabels = storageAxisLabels;
      patch.storageOffset = chunk.offset;
      patch.storageExtent = chunk.extent;

      const size_t axisCount = storageAxisLabels.size();
      std::vector<uint64_t> storageOffsetCanonical(axisCount, 0);
      std::vector<uint64_t> storageExtentCanonical(axisCount, 1);

      for (size_t idx = 0; idx < rawAxisLabels.size(); ++idx) {
        const std::string &label = rawAxisLabels[idx];
        auto iter =
            std::find(storageAxisLabels.begin(), storageAxisLabels.end(), label);
        if (iter == storageAxisLabels.end()) {
          continue;
        }
        size_t storageIndex = static_cast<size_t>(
            std::distance(storageAxisLabels.begin(), iter));
        if (idx < chunk.offset.size()) {
          storageOffsetCanonical[storageIndex] = chunk.offset[idx];
        }
        if (idx < chunk.extent.size()) {
          storageExtentCanonical[storageIndex] = chunk.extent[idx];
        }
      }

      patch.storageOffsetCanonical = storageOffsetCanonical;
      patch.storageExtentCanonical = storageExtentCanonical;

      openPMD::Offset vtkOffset(3, 0);
      openPMD::Extent vtkExtent(3, 1);
      std::array<int, 3> storageToVtk{{-1, -1, -1}};
      for (size_t axis = 0; axis < 3; ++axis) {
        int storageIndex = -1;
        if (axis < storageToVtkVec.size()) {
          storageIndex = storageToVtkVec[axis];
        }
        storageToVtk[axis] = storageIndex;
        if (storageIndex >= 0 &&
            storageIndex < static_cast<int>(storageOffsetCanonical.size())) {
          vtkOffset[axis] =
              storageOffsetCanonical[static_cast<size_t>(storageIndex)];
        }
        if (storageIndex >= 0 &&
            storageIndex < static_cast<int>(storageExtentCanonical.size())) {
          vtkExtent[axis] =
              storageExtentCanonical[static_cast<size_t>(storageIndex)];
        }
      }

      if (chunkCentering == AVT_NODECENT) {
        for (int axis = 0; axis < 3; ++axis) {
          uint64_t localExtent = vtkExtent[axis];
          uint64_t offsetValue = vtkOffset[axis];
          uint64_t globalAxisExtent = globalExtent[axis];
          if (localExtent == 0) {
            localExtent = 1;
          }
          if (offsetValue + localExtent < globalAxisExtent) {
            vtkExtent[axis] = localExtent + 1;
          }
        }
      }

      patch.storageToVtk = storageToVtk;
      patch.offset = vtkOffset;
      patch.extent = vtkExtent;

      for (int axis = 0; axis < 3; ++axis) {
        double spacing = spacingPerAxis[axis];
        patch.spacing[axis] = spacing;

        uint64_t offsetValue = 0;
        if (axis < static_cast<int>(patch.offset.size())) {
          offsetValue = patch.offset[axis];
        }

        double originValue = getValueOr(geom.gridOrigin, axis, 0.0);
        double gridOrigin = unitSI * originValue;
        double positionValue = getValueOr(geom.position, axis, 0.0);

        double logicalIndex = static_cast<double>(offsetValue) + positionValue;
        if (patch.centering == AVT_ZONECENT) {
          logicalIndex -= 0.5;
        }

        patch.origin[axis] = gridOrigin + logicalIndex * spacing;
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

    debug5 << "[openpmd-api-plugin] Nesting domain=" << patchIdx
           << " level=" << levelIndex << " children="
           << JoinContainer(children) << " logicalExtents="
           << JoinContainer(logicalExtents) << '\n';
  }

  debug1 << "[openpmd-api-plugin] BuildDomainNesting completed\n";
  return nesting;
}

avtLocalStructuredDomainBoundaryList *
avtopenpmdFileFormat::BuildDomainBoundaryList(
    const MeshPatchHierarchy &hierarchy, int domain) const {
  if (hierarchy.patches.empty()) {
    return nullptr;
  }

  if (domain < 0 || domain >= static_cast<int>(hierarchy.patches.size())) {
    return nullptr;
  }

  const PatchInfo &patch = hierarchy.patches[domain];
  LogPatchSummary(patch, "BuildDomainBoundaryList patch");

  int extents[6] = {0, 0, 0, 0, 0, 0};
  for (int axis = 0; axis < 3; ++axis) {
    extents[2 * axis] = patch.logicalLower[axis];
    extents[2 * axis + 1] = patch.logicalUpper[axis];
  }

  auto *list = new avtLocalStructuredDomainBoundaryList(domain, extents);

  auto rangesOverlap = [](int a0, int a1, int b0, int b1) {
    return std::max(a0, b0) <= std::min(a1, b1);
  };

  for (size_t otherIdx = 0; otherIdx < hierarchy.patches.size(); ++otherIdx) {
    if (static_cast<int>(otherIdx) == domain) {
      continue;
    }

    const PatchInfo &other = hierarchy.patches[otherIdx];

    int touchAxis = -1;
    int orientation[3] = {0, 0, 0};

    for (int axis = 0; axis < 3; ++axis) {
      bool overlapsOtherDims = true;
      for (int otherAxis = 0; otherAxis < 3; ++otherAxis) {
        if (otherAxis == axis) {
          continue;
        }
        if (!rangesOverlap(patch.logicalLower[otherAxis], patch.logicalUpper[otherAxis],
                           other.logicalLower[otherAxis], other.logicalUpper[otherAxis])) {
          overlapsOtherDims = false;
          break;
        }
      }

      if (!overlapsOtherDims) {
        continue;
      }

      if (patch.logicalUpper[axis] + 1 == other.logicalLower[axis]) {
        if (touchAxis != -1) {
          touchAxis = -2;
          break;
        }
        touchAxis = axis;
        orientation[axis] = 1;
      } else if (other.logicalUpper[axis] + 1 == patch.logicalLower[axis]) {
        if (touchAxis != -1) {
          touchAxis = -2;
          break;
        }
        touchAxis = axis;
        orientation[axis] = -1;
      }
    }

    if (touchAxis == -1) {
      continue;
    }
    if (touchAxis == -2) {
      continue;
    }

    int boundaryExtents[6] = {0, 0, 0, 0, 0, 0};
    for (int axis = 0; axis < 3; ++axis) {
      if (axis == touchAxis) {
        if (orientation[axis] > 0) {
          boundaryExtents[2 * axis] = patch.logicalUpper[axis];
          boundaryExtents[2 * axis + 1] = patch.logicalUpper[axis];
        } else {
          boundaryExtents[2 * axis] = patch.logicalLower[axis];
          boundaryExtents[2 * axis + 1] = patch.logicalLower[axis];
        }
      } else {
        int lower = std::max(patch.logicalLower[axis], other.logicalLower[axis]);
        int upper = std::min(patch.logicalUpper[axis], other.logicalUpper[axis]);
        boundaryExtents[2 * axis] = lower;
        boundaryExtents[2 * axis + 1] = upper;
      }
    }

    list->AddNeighbor(static_cast<int>(otherIdx), static_cast<int>(domain), orientation, boundaryExtents);
  }

  return list;

}

void avtopenpmdFileFormat::AddGhostZonesForPatch(
    const MeshPatchHierarchy &hierarchy, int patchIdx,
    vtkRectilinearGrid *grid) const {
  if (grid == nullptr) {
    return;
  }

  if (patchIdx < 0 || patchIdx >= static_cast<int>(hierarchy.patches.size())) {
    return;
  }

  const PatchInfo &patch = hierarchy.patches[patchIdx];
  if (patch.centering != AVT_ZONECENT) {
    return;
  }

  const int levelIndex = hierarchy.levelIdsPerPatch[patchIdx];
  if (levelIndex + 1 >= hierarchy.numLevels) {
    return;
  }

  const auto &candidateChildren = hierarchy.patchesPerLevel[levelIndex + 1];
  if (candidateChildren.empty()) {
    return;
  }

  std::vector<std::array<int, 6>> refinedRanges;
  refinedRanges.reserve(candidateChildren.size());

  for (int childIdx : candidateChildren) {
    if (childIdx < 0 ||
        childIdx >= static_cast<int>(hierarchy.patches.size())) {
      continue;
    }

    const PatchInfo &child = hierarchy.patches[childIdx];
    if (!IsChildPatch(patch, child)) {
      continue;
    }

    std::array<int, 6> range{0, -1, 0, -1, 0, -1};
    bool hasOverlap = true;

    for (int axis = 0; axis < 3; ++axis) {
      const int parentLower = patch.logicalLower[axis];
      const int parentUpper = patch.logicalUpper[axis];

      const int childLower = child.logicalLower[axis];
      const int childUpper = child.logicalUpper[axis];

      if (parentUpper < parentLower || childUpper < childLower) {
        hasOverlap = false;
        break;
      }

      const double parentSpacing = std::abs(patch.spacing[axis]);
      const double childSpacing = std::abs(child.spacing[axis]);

      int refine = 1;
      if (parentSpacing > 0.0 && childSpacing > 0.0) {
        refine = std::max(1, static_cast<int>(std::lround(parentSpacing / childSpacing)));
      }

      int coarseLower = childLower;
      int coarseUpper = childUpper;
      if (refine > 1) {
        coarseLower = childLower / refine;
        coarseUpper = (childUpper + refine) / refine - 1;
      }

      int overlapLower = std::max(parentLower, coarseLower);
      int overlapUpper = std::min(parentUpper, coarseUpper);
      if (overlapLower > overlapUpper) {
        hasOverlap = false;
        break;
      }

      const int localLower = overlapLower - parentLower;
      const int localUpper = overlapUpper - parentLower;

      range[axis] = localLower;
      range[axis + 3] = localUpper;
    }

    if (hasOverlap) {
      refinedRanges.push_back(range);
    }
  }

  if (refinedRanges.empty()) {
    return;
  }

  const int nx = patch.extent.size() > 0
                     ? static_cast<int>(patch.extent[0])
                     : 1;
  const int ny = patch.extent.size() > 1
                     ? static_cast<int>(patch.extent[1])
                     : 1;
  const int nz = patch.extent.size() > 2
                     ? static_cast<int>(patch.extent[2])
                     : 1;

  const vtkIdType cellCount = static_cast<vtkIdType>(nx) * ny * nz;
  if (cellCount <= 0) {
    return;
  }

  vtkUnsignedCharArray *ghostArray = vtkUnsignedCharArray::New();
  ghostArray->SetNumberOfComponents(1);
  ghostArray->SetNumberOfTuples(cellCount);
  ghostArray->SetName("avtGhostZones");

  unsigned char *values = ghostArray->GetPointer(0);
  std::fill(values, values + cellCount, 0);

  for (const auto &range : refinedRanges) {
    const int x0 = std::max(0, range[0]);
    const int x1 = std::min(nx - 1, range[3]);
    const int y0 = std::max(0, range[1]);
    const int y1 = std::min(ny - 1, range[4]);
    const int z0 = std::max(0, range[2]);
    const int z1 = std::min(nz - 1, range[5]);

    if (x0 > x1 || y0 > y1 || z0 > z1) {
      continue;
    }

    for (int z = z0; z <= z1; ++z) {
      for (int y = y0; y <= y1; ++y) {
        for (int x = x0; x <= x1; ++x) {
          const vtkIdType cellId =
              x + static_cast<vtkIdType>(nx) *
                      (y + static_cast<vtkIdType>(ny) * z);
          avtGhostData::AddGhostZoneType(values[cellId],
                                         REFINED_ZONE_IN_AMR_GRID);
        }
      }
    }
  }

  grid->GetCellData()->AddArray(ghostArray);
  ghostArray->Delete();
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

  uint64_t vtkElementCount = 1;
  for (auto const &nx : patch.extent) {
    vtkElementCount *= (nx == 0 ? 1 : nx);
  }

  const openPMD::Extent &storageExtent =
      patch.storageExtent.empty() ? patch.extent : patch.storageExtent;

  uint64_t storageElementCount = 1;
  for (auto const &nx : storageExtent) {
    storageElementCount *= (nx == 0 ? 1 : nx);
  }

  debug1 << "[openpmd-api-plugin] LoadScalarPatchData elements="
         << vtkElementCount << " storageElements=" << storageElementCount
         << " unitSI=" << rcomp.unitSI()
         << " datatype=" << static_cast<int>(rcomp.getDatatype()) << "\n";

  vtkDataArray *result = nullptr;

  if (rcomp.getDatatype() == openPMD::Datatype::DOUBLE) {
    vtkDoubleArray *data = vtkDoubleArray::New();
    data->SetNumberOfComponents(1);
    data->SetNumberOfTuples(static_cast<vtkIdType>(vtkElementCount));
    double *buffer = data->GetPointer(0);
    try {
      const openPMD::Offset &storageOffset =
          patch.storageOffset.empty() ? patch.offset : patch.storageOffset;
      rcomp.loadChunkRaw(buffer, storageOffset, storageExtent);
      const_cast<openPMD::Series &>(series_).flush();
    } catch (std::exception const &ex) {
      debug1 << "[openpmd-api-plugin] LoadScalarPatchData exception (double): "
             << ex.what() << "\n";
      data->Delete();
      throw;
    }
    ScaleVarData<double>(buffer, storageElementCount, rcomp.unitSI());
    RemapChunkToVTKLayout(buffer, static_cast<size_t>(storageElementCount), patch);
    DuplicateHighEndNodes(patch, buffer);
    result = data;
  } else if (rcomp.getDatatype() == openPMD::Datatype::FLOAT) {
    vtkFloatArray *data = vtkFloatArray::New();
    data->SetNumberOfComponents(1);
    data->SetNumberOfTuples(static_cast<vtkIdType>(vtkElementCount));
    float *buffer = data->GetPointer(0);
    try {
      const openPMD::Offset &storageOffset =
          patch.storageOffset.empty() ? patch.offset : patch.storageOffset;
      rcomp.loadChunkRaw(buffer, storageOffset, storageExtent);
      const_cast<openPMD::Series &>(series_).flush();
    } catch (std::exception const &ex) {
      debug1 << "[openpmd-api-plugin] LoadScalarPatchData exception (float): "
             << ex.what() << "\n";
      data->Delete();
      throw;
    }
    ScaleVarData<float>(buffer, storageElementCount, rcomp.unitSI());
    RemapChunkToVTKLayout(buffer, static_cast<size_t>(storageElementCount), patch);
    DuplicateHighEndNodes(patch, buffer);
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
    openPMD::Iteration const &iteration, const PatchInfo &patch,
    const std::vector<std::string> &components) const {
  debug1 << "[openpmd-api-plugin] LoadVectorPatchData mesh='" << patch.meshName
         << "' components=" << JoinStrings(components) << "' offset="
         << JoinContainer(patch.offset) << " extent="
         << JoinContainer(patch.extent) << "\n";

  if (components.empty()) {
    debug1 << "[openpmd-api-plugin] LoadVectorPatchData invoked with no"
           << " components; returning nullptr\n";
    return nullptr;
  }

  auto meshIt = iteration.meshes.find(patch.meshName);
  if (meshIt == iteration.meshes.end()) {
    debug1 << "[openpmd-api-plugin] LoadVectorPatchData missing mesh '"
           << patch.meshName << "'\n";
    EXCEPTION1(InvalidVariableException, patch.meshName.c_str());
  }

  vtkDataArray *result = nullptr;
  std::vector<vtkDataArray *> loadedComponents;
  loadedComponents.reserve(components.size());

  try {
    for (auto const &componentName : components) {
      vtkDataArray *componentData =
          LoadScalarPatchData(iteration, patch, componentName);
      if (componentData == nullptr) {
        debug1 << "[openpmd-api-plugin] LoadVectorPatchData failed to load"
               << " component '" << componentName << "'\n";
        EXCEPTION1(InvalidVariableException, componentName.c_str());
      }
      loadedComponents.push_back(componentData);
    }

    vtkIdType tupleCount = loadedComponents.front()->GetNumberOfTuples();
    bool hasDoubleComponent = false;
    bool hasFloatComponent = false;

    for (vtkDataArray *component : loadedComponents) {
      if (component->GetNumberOfComponents() != 1) {
        debug1 << "[openpmd-api-plugin] LoadVectorPatchData component has"
               << " unexpected component count="
               << component->GetNumberOfComponents() << "\n";
        EXCEPTION1(InvalidVariableException, patch.meshName.c_str());
      }
      if (component->GetNumberOfTuples() != tupleCount) {
        debug1 << "[openpmd-api-plugin] LoadVectorPatchData tuple count"
               << " mismatch for a component\n";
        EXCEPTION1(InvalidVariableException, patch.meshName.c_str());
      }
      if (vtkDoubleArray::SafeDownCast(component) != nullptr) {
        hasDoubleComponent = true;
      } else if (vtkFloatArray::SafeDownCast(component) != nullptr) {
        hasFloatComponent = true;
      } else {
        debug1 << "[openpmd-api-plugin] LoadVectorPatchData unsupported"
               << " vtkDataArray subtype for component\n";
        EXCEPTION1(InvalidVariableException, patch.meshName.c_str());
      }
    }

    const int numComponents = static_cast<int>(components.size());
    if (hasDoubleComponent) {
      vtkDoubleArray *vectorArray = vtkDoubleArray::New();
      vectorArray->SetNumberOfComponents(numComponents);
      vectorArray->SetNumberOfTuples(tupleCount);
      double *dest = vectorArray->GetPointer(0);

      for (size_t compIdx = 0; compIdx < loadedComponents.size(); ++compIdx) {
        vtkDataArray *component = loadedComponents[compIdx];
        if (auto *doubleComp = vtkDoubleArray::SafeDownCast(component)) {
          const double *src = doubleComp->GetPointer(0);
          for (vtkIdType tuple = 0; tuple < tupleCount; ++tuple) {
            dest[tuple * numComponents + static_cast<int>(compIdx)] = src[tuple];
          }
        } else if (auto *floatComp = vtkFloatArray::SafeDownCast(component)) {
          const float *src = floatComp->GetPointer(0);
          for (vtkIdType tuple = 0; tuple < tupleCount; ++tuple) {
            dest[tuple * numComponents + static_cast<int>(compIdx)] =
                static_cast<double>(src[tuple]);
          }
        }
      }

      result = vectorArray;
    } else if (hasFloatComponent) {
      vtkFloatArray *vectorArray = vtkFloatArray::New();
      vectorArray->SetNumberOfComponents(numComponents);
      vectorArray->SetNumberOfTuples(tupleCount);
      float *dest = vectorArray->GetPointer(0);

      for (size_t compIdx = 0; compIdx < loadedComponents.size(); ++compIdx) {
        auto *floatComp = vtkFloatArray::SafeDownCast(loadedComponents[compIdx]);
        if (floatComp == nullptr) {
          debug1 << "[openpmd-api-plugin] LoadVectorPatchData expected"
                 << " float component but found different type\n";
          EXCEPTION1(InvalidVariableException, patch.meshName.c_str());
        }
        const float *src = floatComp->GetPointer(0);
        for (vtkIdType tuple = 0; tuple < tupleCount; ++tuple) {
          dest[tuple * numComponents + static_cast<int>(compIdx)] = src[tuple];
        }
      }

      result = vectorArray;
    }

    if (result == nullptr) {
      debug1 << "[openpmd-api-plugin] LoadVectorPatchData could not determine"
             << " a suitable output array type\n";
      EXCEPTION1(InvalidVariableException, patch.meshName.c_str());
    }

    for (vtkDataArray *component : loadedComponents) {
      component->Delete();
    }
    loadedComponents.clear();

    debug1 << "[openpmd-api-plugin] LoadVectorPatchData returning array"
           << " tuples=" << tupleCount << " components=" << numComponents
           << "\n";

    return result;
  } catch (...) {
    for (vtkDataArray *component : loadedComponents) {
      component->Delete();
    }
    if (result != nullptr) {
      result->Delete();
    }
    throw;
  }
}

vtkUnstructuredGrid *avtopenpmdFileFormat::EnsureParticleMesh(
    int timeState, const std::string &visitMeshName) {
  if (timeState < 0 ||
      timeState >= static_cast<int>(particleCache_.size())) {
    debug1 << "[openpmd-api-plugin] EnsureParticleMesh invalid timeState="
           << timeState << " cacheSize=" << particleCache_.size() << "\n";
    EXCEPTION1(InvalidVariableException, visitMeshName.c_str());
  }

  auto &speciesCache = particleCache_.at(timeState);
  auto speciesIt = speciesCache.find(visitMeshName);
  if (speciesIt == speciesCache.end()) {
    debug1 << "[openpmd-api-plugin] EnsureParticleMesh missing particle cache"
           << " for mesh '" << visitMeshName << "' timeState=" << timeState
           << "\n";
    EXCEPTION1(InvalidVariableException, visitMeshName.c_str());
  }

  ParticleSpeciesInfo &info = speciesIt->second;

  if (info.meshCached && info.cachedMesh != nullptr) {
    info.cachedMesh->Register(nullptr);
    return info.cachedMesh;
  }

  unsigned long long iterIdx = iterationIndex_.at(timeState);
  debug1 << "[openpmd-api-plugin] EnsureParticleMesh constructing mesh="
         << visitMeshName << " iteration=" << iterIdx << "\n";
  openPMD::Iteration iteration = series_.snapshots()[iterIdx];

  vtkUnstructuredGrid *grid = CreateParticleMesh(iteration, info);
  if (grid == nullptr) {
    debug1 << "[openpmd-api-plugin] CreateParticleMesh returned nullptr for"
           << " mesh '" << visitMeshName << "'\n";
    EXCEPTION1(InvalidVariableException, visitMeshName.c_str());
  }

  info.cachedMesh = grid;
  info.meshCached = true;
  grid->Delete();
  info.cachedMesh->Register(nullptr);
  return info.cachedMesh;
}

vtkUnstructuredGrid *avtopenpmdFileFormat::CreateParticleMesh(
    openPMD::Iteration const &iteration, ParticleSpeciesInfo &info) const {
  auto speciesIt = iteration.particles.find(info.speciesName);
  if (speciesIt == iteration.particles.end()) {
    debug1 << "[openpmd-api-plugin] CreateParticleMesh missing species '"
           << info.speciesName << "'\n";
    return nullptr;
  }

  const openPMD::ParticleSpecies &species = speciesIt->second;
  auto positionIt = species.find("position");
  vtkUnstructuredGrid *grid = vtkUnstructuredGrid::New();

  if (positionIt == species.end()) {
    debug1 << "[openpmd-api-plugin] Species '" << info.speciesName
           << "' has no position record; returning empty point mesh.\n";
    vtkNew<vtkPoints> emptyPoints;
    emptyPoints->SetNumberOfPoints(0);
    grid->SetPoints(emptyPoints);
    return grid;
  }

  const openPMD::Record &positionRecord = positionIt->second;

  std::size_t particleCount = info.particleCount;
  if (particleCount == 0) {
    std::size_t inferred =
        ElementCountFromExtent(positionRecord.getExtent());
    if (inferred == 0 && !info.positionComponents.empty()) {
      const std::string &componentName = info.positionComponents.front();
      if (componentName == openPMD::RecordComponent::SCALAR) {
        inferred = ElementCountFromExtent(positionRecord.getExtent());
      } else {
        auto compIt = positionRecord.find(componentName);
        if (compIt != positionRecord.end()) {
          inferred = ElementCountFromExtent(compIt->second.getExtent());
        }
      }
    }
    particleCount = inferred;
    info.particleCount = particleCount;
  }

  vtkNew<vtkPoints> points;
  points->SetDataTypeToDouble();
  points->SetNumberOfPoints(static_cast<vtkIdType>(particleCount));

  std::array<std::vector<double>, 3> vtkCoords;
  for (auto &axis : vtkCoords) {
    axis.assign(particleCount, 0.0);
  }

  for (std::size_t storageIdx = 0;
       storageIdx < info.positionComponents.size(); ++storageIdx) {
    if (storageIdx >= info.storageToVtk.size()) {
      continue;
    }
    int vtkAxis = info.storageToVtk[storageIdx];
    if (vtkAxis < 0 || vtkAxis >= 3) {
      continue;
    }

    const std::string &componentName = info.positionComponents[storageIdx];
    openPMD::RecordComponent component = positionRecord;
    if (componentName != openPMD::RecordComponent::SCALAR) {
      auto compIt = positionRecord.find(componentName);
      if (compIt == positionRecord.end()) {
        debug1 << "[openpmd-api-plugin] CreateParticleMesh component '"
               << componentName << "' missing for species '"
               << info.speciesName << "'; skipping axis" << vtkAxis << "\n";
        continue;
      }
      component = compIt->second;
    }

    openPMD::Extent extent = component.getExtent();
    if (extent.empty()) {
      extent = openPMD::Extent(1, static_cast<uint64_t>(particleCount));
    }
    openPMD::Offset offset(extent.size(), 0u);
    std::size_t elementCount = ElementCountFromExtent(extent);
    if (elementCount == 0) {
      elementCount = particleCount;
    }

    if (component.getDatatype() == openPMD::Datatype::DOUBLE) {
      std::vector<double> buffer(elementCount, 0.0);
      component.loadChunkRaw(buffer.data(), offset, extent);
      const_cast<openPMD::Series &>(series_).flush();
      ScaleVarData<double>(buffer.data(), elementCount,
                           static_cast<double>(component.unitSI()));
      std::size_t copyCount =
          std::min<std::size_t>(buffer.size(), particleCount);
      std::copy(buffer.begin(), buffer.begin() + copyCount,
                vtkCoords[vtkAxis].begin());
    } else if (component.getDatatype() == openPMD::Datatype::FLOAT) {
      std::vector<float> buffer(elementCount, 0.0f);
      component.loadChunkRaw(buffer.data(), offset, extent);
      const_cast<openPMD::Series &>(series_).flush();
      ScaleVarData<float>(buffer.data(), elementCount,
                          static_cast<float>(component.unitSI()));
      std::size_t copyCount =
          std::min<std::size_t>(buffer.size(), particleCount);
      for (std::size_t idx = 0; idx < copyCount; ++idx) {
        vtkCoords[vtkAxis][idx] = static_cast<double>(buffer[idx]);
      }
    } else {
      debug1 << "[openpmd-api-plugin] CreateParticleMesh unsupported" \
                " datatype="
             << static_cast<int>(component.getDatatype())
             << " for component '" << componentName << "'\n";
    }
  }

  for (vtkIdType idx = 0; idx < static_cast<vtkIdType>(particleCount); ++idx) {
    points->SetPoint(idx, vtkCoords[0][idx], vtkCoords[1][idx],
                     vtkCoords[2][idx]);
  }

  grid->SetPoints(points);
  grid->Allocate(static_cast<vtkIdType>(particleCount));

  for (vtkIdType idx = 0; idx < static_cast<vtkIdType>(particleCount); ++idx) {
    vtkIdType pid = idx;
    grid->InsertNextCell(VTK_VERTEX, 1, &pid);
  }

  debug1 << "[openpmd-api-plugin] CreateParticleMesh built mesh '"
         << info.visitMeshName << "' particles=" << particleCount << "\n";
  return grid;
}

vtkDataArray *avtopenpmdFileFormat::LoadParticleScalarData(
    openPMD::Iteration const &iteration, const ParticleSpeciesInfo &info,
    const std::string &componentPath) const {
  auto speciesIt = iteration.particles.find(info.speciesName);
  if (speciesIt == iteration.particles.end()) {
    debug1 << "[openpmd-api-plugin] LoadParticleScalarData missing species '"
           << info.speciesName << "'\n";
    EXCEPTION1(InvalidVariableException, info.speciesName.c_str());
  }

  const openPMD::ParticleSpecies &species = speciesIt->second;
  auto [recordName, componentName] = SplitRecordComponentPath(componentPath);
  auto recordIt = species.find(recordName);
  if (recordIt == species.end()) {
    debug1 << "[openpmd-api-plugin] LoadParticleScalarData missing record '"
           << recordName << "' for species '" << info.speciesName << "'\n";
    EXCEPTION1(InvalidVariableException, recordName.c_str());
  }

  const openPMD::Record &record = recordIt->second;
  openPMD::RecordComponent rcomp = record;
  if (componentName != openPMD::RecordComponent::SCALAR) {
    auto componentIt = record.find(componentName);
    if (componentIt == record.end()) {
      debug1 << "[openpmd-api-plugin] LoadParticleScalarData missing component '"
             << componentName << "' for record '" << recordName << "'\n";
      EXCEPTION1(InvalidVariableException, componentName.c_str());
    }
    rcomp = componentIt->second;
  }

  openPMD::Extent extent = rcomp.getExtent();
  if (extent.empty()) {
    extent = openPMD::Extent(1, static_cast<uint64_t>(info.particleCount));
  }
  openPMD::Offset offset(extent.size(), 0u);
  std::size_t elementCount = ElementCountFromExtent(extent);
  if (elementCount == 0) {
    elementCount = info.particleCount;
  }

  vtkDataArray *result = nullptr;

  if (rcomp.getDatatype() == openPMD::Datatype::DOUBLE) {
    vtkDoubleArray *data = vtkDoubleArray::New();
    data->SetNumberOfComponents(1);
    data->SetNumberOfTuples(static_cast<vtkIdType>(elementCount));
    double *buffer = data->GetPointer(0);
    rcomp.loadChunkRaw(buffer, offset, extent);
    const_cast<openPMD::Series &>(series_).flush();
    ScaleVarData<double>(buffer, elementCount,
                         static_cast<double>(rcomp.unitSI()));
    result = data;
  } else if (rcomp.getDatatype() == openPMD::Datatype::FLOAT) {
    vtkFloatArray *data = vtkFloatArray::New();
    data->SetNumberOfComponents(1);
    data->SetNumberOfTuples(static_cast<vtkIdType>(elementCount));
    float *buffer = data->GetPointer(0);
    rcomp.loadChunkRaw(buffer, offset, extent);
    const_cast<openPMD::Series &>(series_).flush();
    ScaleVarData<float>(buffer, elementCount,
                        static_cast<float>(rcomp.unitSI()));
    result = data;
  } else {
    debug1 << "[openpmd-api-plugin] LoadParticleScalarData unsupported datatype"
           << "=" << static_cast<int>(rcomp.getDatatype()) << "\n";
  }

  if (result == nullptr) {
    EXCEPTION1(InvalidVariableException, componentPath.c_str());
  }

  debug1 << "[openpmd-api-plugin] LoadParticleScalarData record='" << recordName
         << "' component='" << componentName
         << "' tuples=" << result->GetNumberOfTuples() << "\n";
  return result;
}

vtkDataArray *avtopenpmdFileFormat::LoadParticleVectorData(
    openPMD::Iteration const &iteration, const ParticleSpeciesInfo &info,
    const std::vector<std::string> &componentPaths) const {
  if (componentPaths.empty()) {
    EXCEPTION1(InvalidVariableException, info.speciesName.c_str());
  }

  std::vector<vtkDataArray *> components;
  components.reserve(componentPaths.size());

  try {
    for (auto const &componentPath : componentPaths) {
      vtkDataArray *component =
          LoadParticleScalarData(iteration, info, componentPath);
      components.push_back(component);
    }

    vtkIdType tupleCount = components.front()->GetNumberOfTuples();
    bool hasDouble = false;
    bool hasFloat = false;
    for (vtkDataArray *component : components) {
      if (vtkDoubleArray::SafeDownCast(component) != nullptr) {
        hasDouble = true;
      }
      if (vtkFloatArray::SafeDownCast(component) != nullptr) {
        hasFloat = true;
      }
    }

    vtkDataArray *result = nullptr;
    if (hasDouble || !hasFloat) {
      auto *vec = vtkDoubleArray::New();
      vec->SetNumberOfComponents(static_cast<int>(components.size()));
      vec->SetNumberOfTuples(tupleCount);
      for (std::size_t compIdx = 0; compIdx < components.size(); ++compIdx) {
        vtkDataArray *component = components[compIdx];
        for (vtkIdType tupleIdx = 0; tupleIdx < tupleCount; ++tupleIdx) {
          vec->SetComponent(tupleIdx, static_cast<int>(compIdx),
                            component->GetComponent(tupleIdx, 0));
        }
      }
      result = vec;
    } else {
      auto *vec = vtkFloatArray::New();
      vec->SetNumberOfComponents(static_cast<int>(components.size()));
      vec->SetNumberOfTuples(tupleCount);
      for (std::size_t compIdx = 0; compIdx < components.size(); ++compIdx) {
        vtkDataArray *component = components[compIdx];
        for (vtkIdType tupleIdx = 0; tupleIdx < tupleCount; ++tupleIdx) {
          vec->SetComponent(tupleIdx, static_cast<int>(compIdx),
                            component->GetComponent(tupleIdx, 0));
        }
      }
      result = vec;
    }

    for (vtkDataArray *component : components) {
      component->Delete();
    }

    debug1 << "[openpmd-api-plugin] LoadParticleVectorData components="
           << componentPaths.size() << " tuples=" << tupleCount << "\n";
    return result;
  } catch (...) {
    for (vtkDataArray *component : components) {
      component->Delete();
    }
    throw;
  }
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

  auto meshTypeIt = meshMap_.find(visit_meshname);
  if (meshTypeIt != meshMap_.end() &&
      std::get<0>(meshTypeIt->second) == DatasetType::ParticleSpecies) {
    vtkUnstructuredGrid *grid = EnsureParticleMesh(timeState, visit_meshname);
    if (grid == nullptr) {
      debug1 << "[openpmd-api-plugin] GetMesh failed to build particle mesh '"
             << meshName << "'\n";
      EXCEPTION1(InvalidVariableException, visit_meshname);
    }
    debug1 << "[openpmd-api-plugin] GetMesh success (particle) mesh='"
           << meshName << "'\n";
    return grid;
  }

  EnsureHierarchyInitialized(timeState);

  auto &hierarchyMap = meshHierarchyCache_.at(timeState);
  auto hierarchyMapIt = hierarchyMap.find(visit_meshname);
  if (hierarchyMapIt == hierarchyMap.end()) {
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
  if (auto *rect = vtkRectilinearGrid::SafeDownCast(grid)) {
    AddGhostZonesForPatch(hierarchy, domain, rect);
  }
  debug1 << "[openpmd-api-plugin] GetMesh success mesh='" << meshName
         << "' domain=" << domain << "\n";
  return grid;
}

GeometryData avtopenpmdFileFormat::GetGeometry3D(
    openPMD::Mesh const &mesh,
    openPMD::MeshRecordComponent const *representative,
    bool insertMissingAxes) {
  // get dataOrder
  auto dataOrder = GetEffectiveMeshDataOrder(mesh);
  if (doOverrideMeshDataOrder_) {
    debug5 << "[openpmd-api-plugin] overrideMeshDataOrder applied: "
           << DataOrderToString(dataOrder) << " (original="
           << DataOrderToString(mesh.dataOrder()) << ")\n";
  }

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

  std::vector<double> position;
  bool havePosition = false;
  try {
    position = mesh.position<double>();
    havePosition = !position.empty();
  } catch (openPMD::Error const &) {
    position.clear();
  }
  if (!havePosition && representative != nullptr) {
    try {
      position = representative->position<double>();
      havePosition = !position.empty();
    } catch (openPMD::Error const &) {
      position.clear();
    }
  }
  if (!havePosition) {
    position.assign(axisLabels.size(), 0.0);
  }

  // reverse ordering if mesh.DataOrder() == DataOrder::F
  if (dataOrder == openPMD::Mesh::DataOrder::F) {
    std::reverse(axisLabels.begin(), axisLabels.end());
    std::reverse(extent.begin(), extent.end());
    std::reverse(gridSpacing.begin(), gridSpacing.end());
    std::reverse(gridOrigin.begin(), gridOrigin.end());
    std::reverse(position.begin(), position.end());
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
        position.insert(position.begin(), 0.0);
      }
    }
  }

  GeometryData geom;
  geom.axisLabels = axisLabels;
  geom.storageAxisLabels = axisLabels;
  geom.extent = extent;
  geom.gridOrigin = gridOrigin;
  geom.gridSpacing = gridSpacing;
  if (position.size() != axisLabels.size()) {
    position.resize(axisLabels.size(), 0.0);
  }
  geom.position = position;
  debug2 << "[openpmd-api-plugin] GetGeometry3D result axisLabels="
         << JoinStrings(geom.axisLabels) << " extent="
         << JoinContainer(geom.extent) << " spacing="
         << JoinContainer(geom.gridSpacing) << " origin="
         << JoinContainer(geom.gridOrigin) << " position="
         << JoinContainer(geom.position) << "\n";
  return geom;
}

GeometryData avtopenpmdFileFormat::GetGeometryXYZ(
    openPMD::Mesh const &mesh,
    openPMD::MeshRecordComponent const *representative,
    bool insertMissingAxes) {
  GeometryData geom = GetGeometry3D(mesh, representative, insertMissingAxes);

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

  geom.storageToVtk = transpose;

  // transpose geometry data
  TransposeVector(geom.axisLabels, transpose);
  TransposeVector(geom.extent, transpose);
  TransposeVector(geom.gridSpacing, transpose);
  TransposeVector(geom.gridOrigin, transpose);
  TransposeVector(geom.position, transpose);

  // verify we did it right
  assert(geom.axisLabels[0] == std::string("x"));
  assert(geom.axisLabels[1] == std::string("y"));
  assert(geom.axisLabels[2] == std::string("z"));

  debug2 << "[openpmd-api-plugin] GetGeometryXYZ result axisLabels="
         << JoinStrings(geom.axisLabels) << " extent="
         << JoinContainer(geom.extent) << " spacing="
         << JoinContainer(geom.gridSpacing) << " origin="
         << JoinContainer(geom.gridOrigin) << " position="
         << JoinContainer(geom.position) << "\n";

  return geom;
}

openPMD::Mesh::DataOrder
avtopenpmdFileFormat::GetEffectiveMeshDataOrder(openPMD::Mesh const &mesh) const {
  if (doOverrideMeshDataOrder_) {
    return overrideMeshDataOrder_;
  }
  return mesh.dataOrder();
}

std::vector<int> avtopenpmdFileFormat::GetAxisTranspose(
    std::vector<std::string> const &axisLabelsSrc,
    std::vector<std::string> const &axisLabelsDst) {
  auto getIndexOf = [](std::string const &label,
                       std::vector<std::string> const &axes) {
    auto iter = std::find(axes.begin(), axes.end(), label);
    if (iter != axes.end()) {
      return static_cast<int>(std::distance(axes.begin(), iter));
    }
    return -1;
  };

  std::vector<int> transpose;
  transpose.reserve(axisLabelsDst.size());
  for (auto const &axis : axisLabelsDst) {
    int idx = getIndexOf(axis, axisLabelsSrc);
    if (idx != -1) {
      transpose.push_back(idx);
    }
  }
  return transpose;
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

  auto meshTypeIt = meshMap_.find(visitMeshName);
  if (meshTypeIt != meshMap_.end() &&
      std::get<0>(meshTypeIt->second) == DatasetType::ParticleSpecies) {
    if (timeState < 0 ||
        timeState >= static_cast<int>(particleCache_.size())) {
      debug1 << "[openpmd-api-plugin] GetVar particle cache timeState invalid"
             << " timeState=" << timeState
             << " cacheSize=" << particleCache_.size() << "\n";
      EXCEPTION1(InvalidVariableException, visitMeshName.c_str());
    }

    auto &speciesCache = particleCache_.at(timeState);
    auto speciesIt = speciesCache.find(visitMeshName);
    if (speciesIt == speciesCache.end()) {
      debug1 << "[openpmd-api-plugin] GetVar missing particle metadata for mesh"
             << " '" << visitMeshName << "'\n";
      EXCEPTION1(InvalidVariableException, visitMeshName.c_str());
    }

    unsigned long long iter = iterationIndex_.at(timeState);
    openPMD::Iteration iteration = series_.snapshots()[iter];
    vtkDataArray *data =
        LoadParticleScalarData(iteration, speciesIt->second, component);
    if (data == nullptr) {
      debug1 << "[openpmd-api-plugin] GetVar particle loader returned nullptr"
             << " var='" << requestedVar << "'\n";
      EXCEPTION1(InvalidVariableException, varname);
    }

    debug1 << "[openpmd-api-plugin] GetVar success (particle) var='"
           << requestedVar << "'\n";
    return data;
  }

  EnsureHierarchyInitialized(timeState);

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

vtkDataArray *avtopenpmdFileFormat::GetVectorVar(int timeState, int domain,
                                                 const char *varname) {
  const char *requestedVar = varname != nullptr ? varname : "<null>";
  debug1 << "[openpmd-api-plugin] GetVectorVar timeState=" << timeState
         << " domain=" << domain << " var=" << requestedVar << "\n";

  if (varname == nullptr) {
    debug1 << "[openpmd-api-plugin] GetVectorVar received null var name\n";
    EXCEPTION1(InvalidVariableException, requestedVar);
  }

  auto varIt = vectorVarMap_.find(varname);
  if (varIt == vectorVarMap_.end()) {
    debug1 << "[openpmd-api-plugin] GetVectorVar missing var '" << requestedVar
           << "'\n";
    EXCEPTION1(InvalidVariableException, varname);
  }

  const std::string &visitMeshName = std::get<0>(varIt->second);
  const std::vector<std::string> &components = std::get<1>(varIt->second);

  auto meshTypeIt = meshMap_.find(visitMeshName);
  if (meshTypeIt != meshMap_.end() &&
      std::get<0>(meshTypeIt->second) == DatasetType::ParticleSpecies) {
    if (timeState < 0 ||
        timeState >= static_cast<int>(particleCache_.size())) {
      debug1 << "[openpmd-api-plugin] GetVectorVar particle cache timeState"
             << " invalid timeState=" << timeState
             << " cacheSize=" << particleCache_.size() << "\n";
      EXCEPTION1(InvalidVariableException, visitMeshName.c_str());
    }

    auto &speciesCache = particleCache_.at(timeState);
    auto speciesIt = speciesCache.find(visitMeshName);
    if (speciesIt == speciesCache.end()) {
      debug1 << "[openpmd-api-plugin] GetVectorVar missing particle metadata"
             << " for mesh '" << visitMeshName << "'\n";
      EXCEPTION1(InvalidVariableException, visitMeshName.c_str());
    }

    unsigned long long iter = iterationIndex_.at(timeState);
    openPMD::Iteration iteration = series_.snapshots()[iter];
    vtkDataArray *data = LoadParticleVectorData(iteration, speciesIt->second,
                                                components);
    if (data == nullptr) {
      debug1 << "[openpmd-api-plugin] GetVectorVar particle loader returned"
             << " nullptr var='" << requestedVar << "'\n";
      EXCEPTION1(InvalidVariableException, varname);
    }

    debug1 << "[openpmd-api-plugin] GetVectorVar success (particle) var='"
           << requestedVar << "'\n";
    return data;
  }

  EnsureHierarchyInitialized(timeState);

  auto &hierarchyMap = meshHierarchyCache_.at(timeState);
  auto meshIt = hierarchyMap.find(visitMeshName);
  if (meshIt == hierarchyMap.end()) {
    debug1 << "[openpmd-api-plugin] GetVectorVar missing hierarchy for mesh '"
           << visitMeshName << "'\n";
    EXCEPTION1(InvalidVariableException, visitMeshName.c_str());
  }

  const MeshPatchHierarchy &hierarchy = meshIt->second;
  if (domain < 0 || domain >= static_cast<int>(hierarchy.patches.size())) {
    debug1 << "[openpmd-api-plugin] GetVectorVar invalid domain index "
           << domain << " for mesh '" << visitMeshName << "'\n";
    EXCEPTION1(InvalidVariableException, visitMeshName.c_str());
  }

  const PatchInfo &patch = hierarchy.patches.at(domain);
  std::ostringstream ctx;
  ctx << "GetVectorVar patch domain=" << domain
      << " components=" << JoinStrings(components);
  LogPatchSummary(patch, ctx.str());

  unsigned long long iter = iterationIndex_.at(timeState);
  debug1 << "[openpmd-api-plugin] GetVectorVar loading iteration=" << iter
         << " components='" << JoinStrings(components) << "'\n";
  openPMD::Iteration iteration = series_.snapshots()[iter];

  vtkDataArray *data = LoadVectorPatchData(iteration, patch, components);
  if (data == nullptr) {
    debug1 << "[openpmd-api-plugin] GetVectorVar LoadVectorPatchData returned"
           << " nullptr for var '" << requestedVar << "'\n";
    EXCEPTION1(InvalidVariableException, varname);
  }

  debug1 << "[openpmd-api-plugin] GetVectorVar success var='" << requestedVar
         << "'\n";

  return data;
}

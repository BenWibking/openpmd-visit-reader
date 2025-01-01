// ****************************************************************************
//  dataLayoutTransform.C
// ****************************************************************************

#include <DebugStream.h>
#include <cstddef>
#include <string>

#include "avtopenpmdFileFormat.h"
#include "dataLayoutTransform.h"

std::vector<int> avtopenpmdFileFormat::GetAxisTranspose(
    std::vector<std::string> const &axisLabels) {
  const std::vector<std::string> cartesianAxes = {
      std::string("x"), std::string("y"), std::string("z")};

  auto getIndexOf = [](std::string const &e,
                       std::vector<std::string> const &v) {
    auto iter = std::find(v.begin(), v.end(), e);
    if (iter != v.end()) {
      return std::distance(v.begin(), iter);
    } else {
      return -1L;
    }
  };

  // compute transposition that takes axisLabels to {'x', 'y', 'z'}
  std::vector<int> transpose{};
  for (auto axis : cartesianAxes) {
    auto idx = getIndexOf(axis, axisLabels);
    if (idx != -1L) {
      transpose.push_back(idx);
    }
  }
  return transpose;
}

std::map<std::string, size_t>
GetExtentMap(std::vector<std::string> const &axisLayoutOrder,
             std::vector<uint64_t> const &extents) {
  // return std::map from axis label to extent size
  std::map<std::string, size_t> extent_map{};
  for (int i = 0; i < axisLayoutOrder.size(); ++i) {
    extent_map[axisLayoutOrder[i]] = extents[i];
  }
  return extent_map;
}

std::vector<uint64_t>
GetExtentsWithOrder(std::map<std::string, size_t> &extent_map,
                    std::vector<std::string> const &axisLayoutOrder) {
  // return extents in order
  std::vector<uint64_t> extents(axisLayoutOrder.size());
  for (int i = 0; i < axisLayoutOrder.size(); ++i) {
    extents[i] = extent_map[axisLayoutOrder[i]];
  }
  return extents;
}

std::map<std::string, size_t>
GetLayoutStride(std::vector<std::string> const &axisLayoutOrder,
                std::map<std::string, size_t> &extents) {
  // Returns the strides for a given layout ordering 'axisLayoutOrder'

  // compute coefficients
  // EXAMPLE axisLayoutOrder: 0, 1, 2 == X, Y, Z (from slowest index to fastest)
  // A = 1 * ndims[1] * ndims[2]
  // B = 1 * ndims[2]
  // C = 1;
  // EXAMPLE axisLayoutOrder: 2, 1, 0 == Z, Y, X (from slowest index to fastest)
  // A = 1;
  // B = 1 * ndims[0];
  // C = 1 * ndims[1] * ndims[0];

  // TODO: Is the stride for input and output mixed up??

// [openpmd-api-plugin] GetMesh() for iteration 255 and VisIt mesh rho_mesh
// Mesh transpose: 1, 0, 2
// gridExtent[0] = 51
// gridSpacing[0] = 6e-07
// gridOrigin[0] = -1.5e-05
// gridExtent[1] = 1
// gridSpacing[1] = 0
// gridOrigin[1] = 0
// gridExtent[2] = 201
// gridSpacing[2] = 1e-07
// gridOrigin[2] = 1.02e-05

// [openpmd-api-plugin] Mesh is node-centered.
// [openpmd-api-plugin] GetVar() for iteration 255 and var rho
// Data extents: 51, 201, 
// Data axis labels (slowest-to-fastest layout order): x z 
// Output data extents: 201, 51, 
// VTK axis labels (slowest-to-fastest layout order): z x 
// p = 0
// 	q = 1
// 	stride[x] *= ndims[z]
// p = 1
// p = 0
// 	q = 1
// 	stride[z] *= ndims[x]
// p = 1
// Input layout stride (2D): 1 201
// Output layout stride (2D): 51 1

  std::map<std::string, size_t> stride{};
  for (auto const &axis : axisLayoutOrder) {
    stride[axis] = 1;
  }

  for (int p = 0; p < axisLayoutOrder.size(); ++p) {
    debug5 << "p = " << p << "\n";
    for (int q = p + 1; q < axisLayoutOrder.size(); ++q) {
      debug5 << "\tq = " << q << "\n";
      debug5 << "\tstride[" << axisLayoutOrder[p] << "] *= ndims["
             << axisLayoutOrder[q] << "]\n";

      stride[axisLayoutOrder[p]] *= extents[axisLayoutOrder[q]];
    }
  }

  return stride;
}
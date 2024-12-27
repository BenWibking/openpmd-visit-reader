// ****************************************************************************
//  dataLayoutTransform.C
// ****************************************************************************

#include <cstddef>
#include <cstdint>
#include <string>

#include <DebugStream.h>

#include "avtopenpmdFileFormat.h"

std::vector<int> avtopenpmdFileFormat::GetAxisTranspose(
    std::vector<std::string> const &axisLabels) {
  const std::vector<std::string> cartesianAxes = {
      std::string("x"), std::string("y"), std::string("z")};

  auto getIndexOf = [](std::string const &e,
                       std::vector<std::string> const &v) {
    return std::distance(v.begin(), std::find(v.begin(), v.end(), e));
  };

  // compute transposition that takes axisLabels to {'x', 'y', 'z'}
  std::vector<int> transpose{};
  for (auto axis : cartesianAxes) {
    transpose.push_back(getIndexOf(axis, axisLabels));
  }
  return transpose;
}

std::vector<int>
avtopenpmdFileFormat::GetIndexOrder(std::vector<std::string> axisLabels) {
  // Returns the index ordering (from slowest index to fastest)
  // of data that has data order 'dataOrder' and axis labels 'axisLabels'
  // The index returned for x is 0, y is 1, and z is 2.

  // get index ordering from axisLabels
  std::vector<int> indexOrder{};
  for (std::string const &axisLabel : axisLabels) {
    if (axisLabel == "x") {
      indexOrder.push_back(0);
    } else if (axisLabel == "y") {
      indexOrder.push_back(1);
    } else if (axisLabel == "z") {
      indexOrder.push_back(2);
    }
  }
  return indexOrder;
}

std::tuple<size_t, size_t, size_t>
avtopenpmdFileFormat::GetIndexCoefficients(std::vector<int> const &indexOrder,
                                           std::vector<uint64_t> const &ndims) {
  // Returns the index coefficients A, B, C for a given index ordering

  // compute coefficients
  // EXAMPLE indexOrder: 0, 1, 2 == X, Y, Z (from slowest index to fastest)
  // A = 1 * ndims[1] * ndims[2]
  // B = 1 * ndims[2]
  // C = 1;
  // EXAMPLE indexOrder: 2, 1, 0 == Z, Y, X (from slowest index to fastest)
  // A = 1;
  // B = 1 * ndims[0];
  // C = 1 * ndims[1] * ndims[0];

  // EXAMPLE Data index order: 1 0 2 == Y, X, Z (slowest to fastest)
  // A = 1 * ndims[2]
  // B = 1 * ndims[0] * ndims[2]
  // C = 1
  // EXAMPLE Data index order: 0 2 1  == X, Z, Y (slowest to fastest)
  // A = 1 * ndims[2] * ndims[1] = 201
  // B = 1
  // C = 1 * ndims[1] = 201

  size_t A[3] = {1, 1, 1};
  for (int ii = 0; ii < 3; ++ii) {
    for (int jj = ii + 1; jj < 3; ++jj) {
      A[indexOrder[ii]] *= ndims[indexOrder[jj]];
    }
  }
  return {A[0], A[1], A[2]};
}

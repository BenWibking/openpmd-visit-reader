#ifndef AVT_openpmd_DATA_LAYOUT_TRANSFORM_H
#define AVT_openpmd_DATA_LAYOUT_TRANSFORM_H

#include <DebugStream.h>

#include "avtopenpmdFileFormat.h"

template <typename T>
void avtopenpmdFileFormat::TransposeVector(std::vector<T> &vec_to_transpose,
                                           std::vector<int> const &transpose) {
  // compute transposition that takes axisLabels to {'x', 'y', 'z'}
  std::vector<T> vec = vec_to_transpose;
  for (int i = 0; i < vec.size(); ++i) {
    vec_to_transpose[i] = vec[transpose[i]];
  }
}

template <typename T>
void avtopenpmdFileFormat::TransposeArray(
    T *data_ptr, openPMD::Mesh const &mesh,
    openPMD::Mesh::MeshRecordComponent const &rcomp) {
  /// Transpose data array from input ordering to VTK ordering

  GeometryData geom = GetGeometry3D(mesh);
  auto axisLabels = geom.axisLabels;
  auto extent = geom.extent;

  debug5 << "Data extents: ";
  for (auto const &nx : extent) {
    debug5 << nx << ", ";
  }

  debug5 << "Data axis labels: ";
  for (auto const &label : axisLabels) {
    debug5 << label << " ";
  }
  debug5 << '\n';

  // get index ordering from axisLabels
  std::vector<int> indexOrder = GetIndexOrder(axisLabels);
  debug5 << "Data index order: " << indexOrder[0] << " " << indexOrder[1] << " "
         << indexOrder[2] << '\n';

  // compute the transpose coefficients A, B, C from the axis ordering
  auto [A, B, C] = GetIndexCoefficients(indexOrder, extent);
  debug5 << "Transpose coefficients: A = " << A << ", B = " << B
         << ", C = " << C << '\n';

  // "VTK image data arrays are stored such that the X dimension
  // increases fastest, followed by Y, followed by Z."
  // https://public.kitware.com/pipermail/vtkusers/2016-September/096626.html
  std::vector<int> vtkIndexOrder =
      GetIndexOrder({std::string("z"), std::string("y"), std::string("x")});
  debug5 << "VTK index order: " << vtkIndexOrder[0] << " " << vtkIndexOrder[1]
         << " " << vtkIndexOrder[2] << '\n';

  // transpose coefficients E, F, G
  // VTK indexOrder: 2, 1, 0 == Z, Y, X (from slowest index to fastest)
  // C = ndims[1] * ndims[0] * 1
  // B = ndims[0] * 1
  // A = 1;
  auto [E, F, G] = GetIndexCoefficients(vtkIndexOrder, extent);
  debug5 << "Transpose coefficients: E = " << E << ", F = " << F
         << ", G = " << G << '\n';

  size_t Gp = extent[1] * extent[0] * 1;
  size_t Fp = extent[0] * 1;
  size_t Ep = 1;

  debug5 << "Transpose coefficients: Ep = " << Ep << ", Fp = " << Fp
         << ", Gp = " << Gp << '\n';
  assert(E == Ep);
  assert(F == Fp);
  assert(G == Gp);

  // is this the identity?
  if ((A == E) && (B == F) && (C == G)) {
    // this is the identity transpose, so we don't need to do anything
    debug5 << "TransposeArray(): identity transpose requested, skipping.\n";
    return;
  }

  // do the transpose
  const size_t len = extent[0] * extent[1] * extent[2];
  std::vector<T> data_copy(len);
  std::memcpy(data_copy.data(), data_ptr, len);

  for (int i = 0; i < extent[0]; ++i) {
    for (int j = 0; j < extent[1]; ++j) {
      for (int k = 0; k < extent[2]; ++k) {
        const size_t fromIndex = i * A + j * B + k * C;
        const size_t toIndex = i * E + j * F + k * G;
        assert(fromIndex >= 0);
        assert(fromIndex < len);
        assert(toIndex >= 0);
        assert(toIndex < len);
        data_ptr[toIndex] = data_copy[fromIndex];
      }
    }
  }
}

#endif // AVT_openpmd_DATA_LAYOUT_TRANSFORM_H

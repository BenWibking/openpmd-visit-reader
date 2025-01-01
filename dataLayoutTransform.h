#ifndef AVT_openpmd_DATA_LAYOUT_TRANSFORM_H
#define AVT_openpmd_DATA_LAYOUT_TRANSFORM_H

#include <DebugStream.h>
#include <cstddef>
#include <cstdint>
#include <mdspan/mdspan.hpp>

#include "avtopenpmdFileFormat.h"

std::map<std::string, size_t>
GetExtentMap(std::vector<std::string> const &axisLayoutOrder,
             std::vector<uint64_t> const &extents);

std::vector<uint64_t>
GetExtentsWithOrder(std::map<std::string, size_t> &extent_map,
                    std::vector<std::string> const &axisLayoutOrder);

std::map<std::string, size_t>
GetLayoutStride(std::vector<std::string> const &axisLayoutOrder,
                std::map<std::string, size_t> &extents);

template <typename T>
void avtopenpmdFileFormat::TransposeVector(std::vector<T> &vec_to_transpose,
                                           std::vector<int> const &transpose) {
  // compute transposition of 'vec_to_transpose' according to 'transpose'
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

  GeometryData geom = GetGeometry3D(mesh, false);
  auto axisLabels = geom.axisLabels;
  auto input_extent = geom.extent;

  // get extent map
  auto extent_map = GetExtentMap(axisLabels, input_extent);

  // "VTK image data arrays are stored such that the X dimension increases
  // fastest, followed by Y, followed by Z."
  std::vector<std::string> vtkLabelOrderSlowToFast = {"z", "x"};

  // get output extents (ordered following vtkLabelOrderSlowToFast)
  std::vector<uint64_t> output_extent =
      GetExtentsWithOrder(extent_map, vtkLabelOrderSlowToFast);

  // debugging
  {
    debug5 << "Data extents: ";
    for (auto const &nx : input_extent) {
      debug5 << nx << ", ";
    }
    debug5 << '\n';

    debug5 << "Data axis labels (slowest-to-fastest layout order): ";
    for (auto const &label : axisLabels) {
      debug5 << label << " ";
    }
    debug5 << '\n';

    debug5 << "Output data extents: ";
    for (auto const &nx : output_extent) {
      debug5 << nx << ", ";
    }
    debug5 << '\n';

    debug5 << "VTK axis labels (slowest-to-fastest layout order): ";
    for (auto const &label : vtkLabelOrderSlowToFast) {
      debug5 << label << " ";
    }
    debug5 << '\n';
  }

  // [openpmd-api-plugin] GetMesh() for iteration 255 and VisIt mesh rho_mesh
  // gridExtent[0] = 51
  // gridSpacing[0] = 6e-07
  // gridOrigin[0] = -1.5e-05
  // gridExtent[1] = 201
  // gridSpacing[1] = 1e-07
  // gridOrigin[1] = 1.02e-05

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

  // copy input array
  size_t len = 1;
  for (auto const &nx : input_extent) {
    len *= nx;
  }
  std::vector<T> data_copy(len);
  std::memcpy(data_copy.data(), data_ptr, len);
  std::memset(data_ptr, 0, len); // set to zero

  if (input_extent.size() == 3) {
    /// 3D array transpose

    // extents
    Kokkos::dextents<std::size_t, 3> extent{input_extent[0], input_extent[1],
                                            input_extent[2]};

    // input stride
    auto stride_map = GetLayoutStride(axisLabels, extent_map);
    std::array<size_t, 3> in_stride{stride_map["z"], stride_map["y"],
                                    stride_map["x"]};

    // output stride
    auto vtkStride_map = GetLayoutStride(vtkLabelOrderSlowToFast, extent_map);
    std::array<size_t, 3> out_stride{vtkStride_map["z"], vtkStride_map["y"],
                                     vtkStride_map["x"]};

    debug5 << "Input layout stride (3D): " << in_stride[0] << " "
           << in_stride[1] << " " << in_stride[2] << '\n';
    debug5 << "Output layout stride (3D): " << out_stride[0] << " "
           << out_stride[1] << " " << out_stride[2] << '\n';

    Kokkos::mdspan input_md{data_copy.data(),
                            Kokkos::layout_stride::mapping{extent, in_stride}};
    Kokkos::mdspan output_md{
        data_ptr, Kokkos::layout_stride::mapping{extent, out_stride}};

    for (std::size_t i = 0; i != output_md.extent(0); i++) {
      for (std::size_t j = 0; j != output_md.extent(1); j++) {
        for (std::size_t k = 0; k != output_md.extent(2); k++) {
          // bounds-checking with .at()
          output_md.at(i, j, k) = input_md.at(i, j, k);
        }
      }
    }
  } else if (input_extent.size() == 2) {
    /// 2D array transpose [currently hard-coded to x-z axes]

    // extents
    Kokkos::dextents<std::size_t, 2> extent{input_extent[0], input_extent[1]};

    // input stride
    auto stride_map = GetLayoutStride(axisLabels, extent_map);
    std::array<size_t, 2> in_stride{stride_map["z"], stride_map["x"]};

    // output stride
    auto vtkStride_map = GetLayoutStride(vtkLabelOrderSlowToFast, extent_map);
    std::array<size_t, 2> out_stride{vtkStride_map["z"], vtkStride_map["x"]};

    debug5 << "Input layout stride (2D): " << in_stride[0] << " "
           << in_stride[1] << '\n';
    debug5 << "Output layout stride (2D): " << out_stride[0] << " "
           << out_stride[1] << '\n';

    Kokkos::mdspan input_md{data_copy.data(),
                            Kokkos::layout_stride::mapping{extent, in_stride}};
    Kokkos::mdspan output_md{
        data_ptr, Kokkos::layout_stride::mapping{extent, out_stride}};

    for (std::size_t i = 0; i != output_md.extent(0); i++) {
      for (std::size_t j = 0; j != output_md.extent(1); j++) {
        // bounds-checking with .at()
        output_md.at(i, j) = input_md.at(i, j);
      }
    }
  }
}

#endif // AVT_openpmd_DATA_LAYOUT_TRANSFORM_H

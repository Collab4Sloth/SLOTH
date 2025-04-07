/**
 * @file main.cpp
 * @author cp273896 (clement.plumecocq@cea.fr)
 * @brief
 * @version 0.1
 * @date 2025-04-04
 *
 * @copyright Copyright (c) 2025
 *
 */
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "kernel/sloth.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
#include "tests/tests.hpp"

///---------------
/// Main program
///---------------
int main(int argc, char* argv[]) {
  //---------------------------------------
  // Initialize MPI and HYPRE
  //---------------------------------------

  mfem::Mpi::Init(argc, argv);
  mfem::Hypre::Init();

  int rank = mfem::Mpi::WorldRank();
  if (rank == 0) {
    // 1D verification
    {
      std::array<double, 1> point_to_interpolate = {0.5};
      std::array<std::vector<double>, 1> grid_values = {
          std::vector<double>{-1.0, 0.0, 1., 1.5, 2.0},
      };
      std::array<std::size_t, 1> lower_indices = {1};
      boost::multi_array<double, 1> array(boost::extents[5]);
      array[0] = 1.0;
      array[1] = 2.0;
      array[2] = 3.0;
      array[3] = 4.0;
      array[4] = 5.0;

      std::array<double, 1> alpha;
      alpha = MultiLinearInterpolator<double, 1>::computeInterpolationCoefficients(point_to_interpolate, lower_indices,
                                                             grid_values);

      double result = MultiLinearInterpolator<double, 1>::computeInterpolation(lower_indices, alpha, array,[](double v) { return v; });
      int i = lower_indices[0];
      double u = (point_to_interpolate[0] - grid_values[0][i]) /
                 (grid_values[0][i + 1] - grid_values[0][i]);
      double ref = array[i] + u * (array[i + 1] - array[i]);
      bool isConsistent = std::abs(ref - result) < 1e-15;
      MFEM_VERIFY(isConsistent, "Linear interpolation is not consistent")
    }
    {  // 2D verification
      std::array<double, 2> point_to_interpolate2D = {0.5, 0.5};
      std::array<std::size_t, 2> lower_indices2D = {0, 0};
      std::array<std::vector<double>, 2> grid_values2D = {std::vector<double>{0.0, 1.},
                                                          std::vector<double>{0.0, 1.}};
      boost::multi_array<double, 2> array2D(boost::extents[2][2]);
      array2D[0][0] = 1.;
      array2D[0][1] = 3.;
      array2D[1][0] = 2.;
      array2D[1][1] = 4.;
      std::array<double, 2> alpha2D;
      alpha2D = MultiLinearInterpolator<double, 2>::computeInterpolationCoefficients(point_to_interpolate2D,
                                                                 lower_indices2D, grid_values2D);
      double result = MultiLinearInterpolator<double, 2>::computeInterpolation(lower_indices2D, alpha2D, array2D,[](double v) { return v; });

      int i = lower_indices2D[0];
      int j = lower_indices2D[1];
      double t = (point_to_interpolate2D[0] - grid_values2D[0][i]) /
                 (grid_values2D[0][i + 1] - grid_values2D[0][i]);
      double u = (point_to_interpolate2D[1] - grid_values2D[1][j]) /
                 (grid_values2D[1][j + 1] - grid_values2D[1][j]);
      double ref = (1 - t) * (1 - u) * array2D[i][j] + t * (1 - u) * array2D[i + 1][j] +
                   t * u * array2D[i + 1][j + 1] + (1 - t) * u * array2D[i][j + 1];
      bool isConsistent = std::abs(ref - result) < 1e-15;
      MFEM_VERIFY(isConsistent, "Bilinear interpolation is not consistent")
    }
    {
      // 3D verification
      std::array<double, 3> point_to_interpolate = {0.5, 0.6, 0.7};
      std::array<std::size_t, 3> lower_indices = {0, 0, 0};
      std::array<std::vector<double>, 3> grid_values = {std::vector<double>{0.0, 1.0},
                                                        std::vector<double>{0.0, 1.0},
                                                        std::vector<double>{0.0, 1.0}};

      boost::multi_array<double, 3> array(boost::extents[2][2][2]);
      array[0][0][0] = 1.0;
      array[1][0][0] = 2.0;
      array[0][1][0] = 3.0;
      array[1][1][0] = 4.0;
      array[0][0][1] = 5.0;
      array[1][0][1] = 6.0;
      array[0][1][1] = 7.0;
      array[1][1][1] = 8.0;
      std::array<double, 3> alpha;
      alpha = MultiLinearInterpolator<double, 3>::computeInterpolationCoefficients(point_to_interpolate, lower_indices,
                                                             grid_values);

      double result = MultiLinearInterpolator<double, 3>::computeInterpolation(lower_indices, alpha, array,[](double v) { return v; });

      int i = lower_indices[0];
      int j = lower_indices[1];
      int k = lower_indices[2];
      double t = (point_to_interpolate[0] - grid_values[0][i]) /
                 (grid_values[0][i + 1] - grid_values[0][i]);
      double u = (point_to_interpolate[1] - grid_values[1][j]) /
                 (grid_values[1][j + 1] - grid_values[1][j]);
      double v = (point_to_interpolate[2] - grid_values[2][k]) /
                 (grid_values[2][k + 1] - grid_values[2][k]);

      double ref =
          (1 - t) * (1 - u) * (1 - v) * array[i][j][k] +
          t * (1 - u) * (1 - v) * array[i + 1][j][k] + (1 - t) * u * (1 - v) * array[i][j + 1][k] +
          (1 - t) * (1 - u) * v * array[i][j][k + 1] + t * u * (1 - v) * array[i + 1][j + 1][k] +
          t * (1 - u) * v * array[i + 1][j][k + 1] + (1 - t) * u * v * array[i][j + 1][k + 1] +
          t * u * v * array[i + 1][j + 1][k + 1];
      bool isConsistent = std::abs(ref - result) < 1e-15;
      MFEM_VERIFY(isConsistent, "Bilinear interpolation is not consistent")
    }
  }

  MPI_Finalize();

  return 0;
}

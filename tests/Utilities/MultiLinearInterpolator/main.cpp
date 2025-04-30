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
  // 1D verification
  {
    std::size_t dim = 1;
    std::vector<double> point_to_interpolate = {0.5};
    std::vector<std::vector<double>> grid_values = {
        std::vector<double>{-1.0, 0.0, 1., 1.5, 2.0},
    };
    std::vector<std::size_t> lower_indices = {1};

    std::vector<double> data(5);
    data[0] = 1.0;
    data[1] = 2.0;
    data[2] = 3.0;
    data[3] = 4.0;
    data[4] = 5.0;
    FlattenedTensor<double> farray({1., 2., 3., 4., 5});
    farray.set_dim(1);
    std::vector<std::size_t> shape = {5};
    farray.set_shape(shape);
    std::vector<double> alpha(dim);
    alpha = MultiLinearInterpolator<double>::computeInterpolationCoefficients(
        dim, point_to_interpolate, lower_indices, grid_values);
    double result =
        MultiLinearInterpolator<double>::computeInterpolation(dim, lower_indices, alpha, farray);

    int i = lower_indices[0];
    double u =
        (point_to_interpolate[0] - grid_values[0][i]) / (grid_values[0][i + 1] - grid_values[0][i]);
    double ref = farray[i] + u * (farray[i + 1] - farray[i]);

    bool isConsistent = std::abs(ref - result) < 1e-15;
    MFEM_VERIFY(isConsistent, "Linear interpolation is not consistent")
  }
  {  // 2D verification
    std::size_t dim = 2;
    std::vector<double> point_to_interpolate2D = {0.5, 0.5};
    std::vector<std::size_t> lower_indices2D = {0, 0};
    std::vector<std::vector<double>> grid_values2D = {std::vector<double>{0.0, 1.},
                                                      std::vector<double>{0.0, 1.}};

    FlattenedTensor<double> farray2D({1., 3., 2., 4.});
    farray2D.set_dim(2);
    std::vector<std::size_t> shape = {2, 2};
    farray2D.set_shape(shape);

    std::vector<double> alpha2D(2);
    alpha2D = MultiLinearInterpolator<double>::computeInterpolationCoefficients(
        dim, point_to_interpolate2D, lower_indices2D, grid_values2D);
    double result = MultiLinearInterpolator<double>::computeInterpolation(dim, lower_indices2D,
                                                                          alpha2D, farray2D);

    std::size_t i = lower_indices2D[0];
    std::size_t j = lower_indices2D[1];
    double t = (point_to_interpolate2D[0] - grid_values2D[0][i]) /
               (grid_values2D[0][i + 1] - grid_values2D[0][i]);
    double u = (point_to_interpolate2D[1] - grid_values2D[1][j]) /
               (grid_values2D[1][j + 1] - grid_values2D[1][j]);
    double ref = (1 - t) * (1 - u) * farray2D.evaluate({i, j}) +
                 t * (1 - u) * farray2D.evaluate({i + 1, j}) +
                 t * u * farray2D.evaluate({i + 1, j + 1}) +
                 (1 - t) * u * farray2D.evaluate({i, j + 1});
    bool isConsistent = std::abs(ref - result) < 1e-15;
    MFEM_VERIFY(isConsistent, "Bilinear interpolation is not consistent")
  }
  {
    // 3D verification
    std::size_t dim = 3;
    std::vector<double> point_to_interpolate = {0.5, 0.6, 0.7};
    std::vector<std::size_t> lower_indices = {0, 0, 0};
    std::vector<std::vector<double>> grid_values = {std::vector<double>{0.0, 1.0},
                                                    std::vector<double>{0.0, 1.0},
                                                    std::vector<double>{0.0, 1.0}};


    FlattenedTensor<double> farray;
    farray.resize(8);
    farray.set_dim(3);
    std::vector<std::size_t> shape = {2, 2, 2};
    farray.set_shape(shape);

    farray[farray.flattened_index({0, 0, 0})] = 1.0;
    farray[farray.flattened_index({1, 0, 0})] = 2.0;
    farray[farray.flattened_index({0, 1, 0})] = 3.0;
    farray[farray.flattened_index({1, 1, 0})] = 4.0;
    farray[farray.flattened_index({0, 0, 1})] = 5.0;
    farray[farray.flattened_index({1, 0, 1})] = 6.0;
    farray[farray.flattened_index({0, 1, 1})] = 7.0;
    farray[farray.flattened_index({1, 1, 1})] = 8.0;

    std::vector<double> alpha(dim);
    alpha = MultiLinearInterpolator<double>::computeInterpolationCoefficients(
        dim, point_to_interpolate, lower_indices, grid_values);

    double result =
        MultiLinearInterpolator<double>::computeInterpolation(dim, lower_indices, alpha, farray);

    std::size_t i = lower_indices[0];
    std::size_t j = lower_indices[1];
    std::size_t k = lower_indices[2];
    double t =
        (point_to_interpolate[0] - grid_values[0][i]) / (grid_values[0][i + 1] - grid_values[0][i]);
    double u =
        (point_to_interpolate[1] - grid_values[1][j]) / (grid_values[1][j + 1] - grid_values[1][j]);
    double v =
        (point_to_interpolate[2] - grid_values[2][k]) / (grid_values[2][k + 1] - grid_values[2][k]);

    double ref = (1 - t) * (1 - u) * (1 - v) * farray.evaluate({i, j, k}) +
                 t * (1 - u) * (1 - v) * farray.evaluate({i + 1, j, k}) +
                 (1 - t) * u * (1 - v) * farray.evaluate({i, j + 1, k}) +
                 (1 - t) * (1 - u) * v * farray.evaluate({i, j, k + 1}) +
                 t * u * (1 - v) * farray.evaluate({i + 1, j + 1, k}) +
                 t * (1 - u) * v * farray.evaluate({i + 1, j, k + 1}) +
                 (1 - t) * u * v * farray.evaluate({i, j + 1, k + 1}) +
                 t * u * v * farray.evaluate({i + 1, j + 1, k + 1});
    bool isConsistent = std::abs(ref - result) < 1e-15;
    MFEM_VERIFY(isConsistent, "Bilinear interpolation is not consistent")
  }

  MPI_Finalize();

  return 0;
}

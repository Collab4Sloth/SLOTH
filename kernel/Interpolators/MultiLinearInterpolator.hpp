
/**
 * @file MultiLinearInterpolator.cpp
 * @author cp273896  (clement.plumecocq@cea.fr)
 * @brief Class for defining a multilinear interpolation operator for an abstract dimension
 * @version 0.1
 * @date 2025-04-02
 *
 * Copyright CEA (c) 2025
 *
 */
#include <boost/multi_array.hpp>
#include <functional>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "kernel/Utils/Utils.hpp"

#pragma once
/**
 * @brief 
 * 
 * @tparam T 
 */
template <typename T>
class MultiLinearInterpolator {
 private:
 public:
  MultiLinearInterpolator();
  static double computeInterpolation(const std::size_t &N,
                                     const std::vector<std::size_t> &lower_indices,
                                     const std::vector<double> &alpha,
                                     const FlattenedTensor<double> &array);
  static std::vector<double> computeInterpolationCoefficients(
      const std::size_t &N, const std::vector<double> &point_to_interpolate,
      const std::vector<std::size_t> &lower_indices,
      const std::vector<std::vector<double>> &grid_values);
  ~MultiLinearInterpolator();
};
template <typename T>
MultiLinearInterpolator<T>::MultiLinearInterpolator() {}

/**
 * @brief 
 * 
 * @tparam T 
 * @param N 
 * @param lower_indices 
 * @param alpha 
 * @param array 
 * @return double 
 */
template <typename T>
double MultiLinearInterpolator<T>::computeInterpolation(
    const std::size_t &N, const std::vector<std::size_t> &lower_indices,
    const std::vector<double> &alpha, const FlattenedTensor<double> &array) {
  double val_interpolated = 0.0;
  std::size_t number_of_hypercube_neighbors(1 << N);
  std::size_t vectorized_ind = 0;
  std::vector<std::size_t> indices(N);

  for (std::size_t i = 0; i < number_of_hypercube_neighbors; ++i) {
    double weight = 1.0;

    for (std::size_t j = 0; j < N; ++j) {
      bool condition = (i & (1 << j)) != 0;
      weight *= condition ? alpha[j] : (1.0 - alpha[j]);
      indices[j] = lower_indices[j] + condition;
    }

    // vectorized_ind = array.flattened_index(indices);
    val_interpolated += weight * array.evaluate(indices);
  }
  return val_interpolated;
}

/**
 * @brief 
 * 
 * @tparam T 
 * @param N 
 * @param point_to_interpolate 
 * @param lower_indices 
 * @param grid_values 
 * @return std::vector<double> 
 */
template <typename T>
std::vector<double> MultiLinearInterpolator<T>::computeInterpolationCoefficients(
    const std::size_t &N, const std::vector<double> &point_to_interpolate,
    const std::vector<std::size_t> &lower_indices,
    const std::vector<std::vector<double>> &grid_values) {
  double val_interpolated = 0.0;

  std::vector<double> alpha(N);
  for (std::size_t i = 0; i < N; ++i) {
    double x0 = grid_values[i][lower_indices[i]];
    double x1 = grid_values[i][lower_indices[i] + 1];
    alpha[i] = (point_to_interpolate[i] - x0) / (x1 - x0);
  }
  return alpha;
}

/**
 * @brief Destroy the Multi Linear Interpolator< T>:: Multi Linear Interpolator object
 * 
 * @tparam T 
 */
template <typename T>
MultiLinearInterpolator<T>::~MultiLinearInterpolator() {}

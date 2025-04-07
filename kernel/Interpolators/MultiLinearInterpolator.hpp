
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
#include <functional>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <boost/multi_array.hpp>

#pragma once

template <typename T, std::size_t N>
class MultiLinearInterpolator {
 private:
 public:
  MultiLinearInterpolator();
  template <typename SCALEF>
  static double computeInterpolation(const std::array<std::size_t, N> &lower_indices,
                                     const std::array<double, N> &alpha,
                                     const boost::multi_array<T, N> &array, SCALEF scalling_func);
  static std::array<double, N> computeInterpolationCoefficients(
      const std::array<double, N> &point_to_interpolate,
      const std::array<std::size_t, N> &lower_indices,
      const std::array<std::vector<double>, N> &grid_values);
  ~MultiLinearInterpolator();
};
template <typename T, std::size_t N>
MultiLinearInterpolator<T, N>::MultiLinearInterpolator() {}

/**
 * @brief Calculate interpolation using the alpha weighting factor
 *
 * @tparam T
 * @tparam N
 * @param alpha
 * @param array
 * @param scalling_func
 */
template <typename T, std::size_t N>
template <typename SCALEF>
double MultiLinearInterpolator<T, N>::computeInterpolation(
    const std::array<std::size_t, N> &lower_indices, const std::array<double, N> &alpha,
    const boost::multi_array<T, N> &array, SCALEF scalling_func) {
  double val_interpolated = 0.0;
  std::size_t number_of_hypercube_neighbors(1 << N);

  for (std::size_t i = 0; i < number_of_hypercube_neighbors; ++i) {
    std::array<std::size_t, N> indices;
    double weight = 1.0;

    for (std::size_t j = 0; j < N; ++j) {
      bool condition = (i & (1 << j)) != 0;
      weight *= condition ? alpha[j] : (1.0 - alpha[j]);
      indices[j] = lower_indices[j] + condition;
    }
    val_interpolated += weight * scalling_func(array(indices));
  }

  return val_interpolated;
}
/**
 * @brief Calculation of the alpha weighting factor for interpolation
 *
 * @tparam T
 * @tparam N
 * @param point_to_interpolate
 * @param lower_indices
 * @param grid_values
 */
template <typename T, std::size_t N>
std::array<double, N> MultiLinearInterpolator<T, N>::computeInterpolationCoefficients(
    const std::array<double, N> &point_to_interpolate,
    const std::array<std::size_t, N> &lower_indices,
    const std::array<std::vector<double>, N> &grid_values) {
  double val_interpolated = 0.0;

  std::array<double, N> alpha;
  for (std::size_t i = 0; i < N; ++i) {
    double x0 = grid_values[i][lower_indices[i]];
    double x1 = grid_values[i][lower_indices[i] + 1];
    alpha[i] = (point_to_interpolate[i] - x0) / (x1 - x0);
  }
  return alpha;
}

template <typename T, std::size_t N>
MultiLinearInterpolator<T, N>::~MultiLinearInterpolator() {}

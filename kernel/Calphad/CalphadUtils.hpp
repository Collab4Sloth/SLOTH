/**
 * @file CalphadUtils.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Class of methods usefull for to compute equilibrium state
 * object
 * @version 0.1
 * @date 2025-09-05
 * 
 * Copyright CEA (C) 2025
 * 
 * This file is part of SLOTH.
 * 
 * SLOTH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * SLOTH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

#include <algorithm>
#include <iterator>
#include <string>
#include <tuple>
#include <vector>
#pragma once

template <typename T>
class CalphadUtils {
 public:
  CalphadUtils();

  std::vector<int> sort_nodes(const T &temp, const T &press,
                              const std::string &temperature_sort_method,
                              const std::string &pressure_sort_method);

  size_t get_size(const T &v);

  ~CalphadUtils();
};

/**
 * @brief Construct a new CalphadUtils::CalphadUtils object
 *
 */
template <typename T>
CalphadUtils<T>::CalphadUtils() {}

/**
 * @brief Return the size of the T element
 *
 * @tparam T
 * @param v
 * @return size_t
 */
template <typename T>
size_t CalphadUtils<T>::get_size(const T &v) {
  size_t vsize = 0;
  for (const auto &vel : v) {
    ++vsize;
  }
  return vsize;
}

/**
 * @brief Method used to sort nodes as function of temperature and/or pressure
 *
 * @tparam T
 * @param temp
 * @param press
 * @param temperature_sort_method
 * @param pressure_sort_method
 * @return std::vector<int>
 */
template <typename T>
std::vector<int> CalphadUtils<T>::sort_nodes(const T &temp, const T &press,
                                             const std::string &temperature_sort_method,
                                             const std::string &pressure_sort_method) {
  // ===============
  //  Initialization
  // ===============
  std::vector<std::tuple<int, double, double>> sorted_data;

  size_t nb_nodes = this->get_size(temp);

  for (auto i = 0; i < nb_nodes; i++) {
    sorted_data.emplace_back(std::make_tuple(i, temp[i], press[i]));
  }
  // ===============
  //  Pressure
  // ===============
  switch (pressure_sort_method::from(pressure_sort_method)) {
    case pressure_sort_method::descending: {
      std::sort(
          begin(sorted_data), end(sorted_data),
          [](const std::tuple<int, double, double> &t1, const std::tuple<int, double, double> &t2) {
            return std::get<2>(t1) > std::get<2>(t2);
          });
      break;
    }
    case pressure_sort_method::ascending: {
      std::sort(
          begin(sorted_data), end(sorted_data),
          [](const std::tuple<int, double, double> &t1, const std::tuple<int, double, double> &t2) {
            return std::get<2>(t1) < std::get<2>(t2);
          });
      break;
    }
    case pressure_sort_method::no: {
      // Nothing to done
      break;
    }
    default: {
      throw std::runtime_error(
          "CalphadUtils::sort_nodes: error in the choice of method used to sort "
          "nodes in function of pressure. Available choices are : ascending,  "
          "descending "
          "(default), No");
    }
  }

  // ===============
  //  Temperature
  // ===============
  switch (temperature_sort_method::from(temperature_sort_method)) {
    case temperature_sort_method::descending: {
      std::stable_sort(
          begin(sorted_data), end(sorted_data),
          [](const std::tuple<int, double, double> &t1, const std::tuple<int, double, double> &t2) {
            return std::get<1>(t1) > std::get<1>(t2);
          });
      break;
    }
    case temperature_sort_method::ascending: {
      std::stable_sort(
          begin(sorted_data), end(sorted_data),
          [](const std::tuple<int, double, double> &t1, const std::tuple<int, double, double> &t2) {
            return std::get<1>(t1) < std::get<1>(t2);
          });
      break;
    }
    case temperature_sort_method::no: {
      // Nothing to done
      break;
    }
    default: {
      throw std::runtime_error(
          "CalphadUtils::sort_nodes: error in the choice of method used to sort "
          "nodes in function of temperature. Available choices are : ascending, "
          "descending (default), No");
    }
  }
  std::vector<int> sorted_nodes;
  for (const auto &i : sorted_data) {
    sorted_nodes.emplace_back(std::get<0>(i));
  }
  sorted_data.clear();

  return sorted_nodes;
}

/**
 * @brief Destroy the CalphadUtils::CalphadUtils object
 *
 */
template <typename T>
CalphadUtils<T>::~CalphadUtils() {}

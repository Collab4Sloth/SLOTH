/**
 * @file FittedLoading.hpp
 * @author cp273896 (clement.plumecocq@cea.fr)
 * @brief
 * @version 0.1
 * @date 2025-11-28
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
 *
 */
#include <algorithm>
#include <functional>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "ExternalLoading/LoadingOptions.hpp"
#include "Parameters/Parameter.hpp"

#pragma once

template <CoordinateSystem CS>
class FittedLoading {
 public:
  std::vector<std::function<double(const std::vector<double>&, double)>> funcs;

  const Parameters& params_;
  void do_time_step(double& time, std::vector<mfem::Vector>& unknown_,
                    const std::vector<std::vector<std::string>>& unks_info,
                    const std::vector<std::tuple<std::string, mfem::Vector>>& coordinates_gf);
  void initialize();
  void get_parameters();

  FittedLoading(const Parameters& params);
};

/**
 * @brief Construct a new Fitted Loading<CS>:: Fitted Loading object
 *
 * @tparam CS
 * @param params
 */
template <CoordinateSystem CS>
FittedLoading<CS>::FittedLoading(const Parameters& params) : params_(params) {
  this->get_parameters();
}

/**
 * @brief
 *
 * @tparam CS
 */
template <CoordinateSystem CS>
void FittedLoading<CS>::get_parameters() {
  this->funcs = this->params_.template get_param_value_or_default<
      std::vector<std::function<double(const std::vector<double>&, double)>>>("LoadingFunctions",
                                                                              {});
}

/**
 * @brief
 *
 * @tparam CS
 * @param time
 * @param unknown_
 * @param unks_info
 * @param coordinates_gf
 */
template <CoordinateSystem CS>
void FittedLoading<CS>::do_time_step(
    double& time, std::vector<mfem::Vector>& unknown_,
    const std::vector<std::vector<std::string>>& unks_info,
    const std::vector<std::tuple<std::string, mfem::Vector>>& coordinates_gf) {
  std::vector<double> coord;

  for (std::size_t u = 0; u < unknown_.size(); u++) {  // loop on unknown
    auto& unk = unknown_[u];
    for (std::size_t i = 0; i < unk.Size(); i++) {  // loop on each point of the unknows

      std::size_t dim = coordinates_gf.size();  // problem dimension

      if constexpr (CS == CoordinateSystem::Cartesian) {
        coord.resize(dim);
        for (std::size_t j = 0; j < dim; j++) {
          coord[j] = std::get<1>(coordinates_gf[j])(i);  // get physical coordinates (x,y,z)
        }
      } else if constexpr (CS == CoordinateSystem::Axisymmetric) {
        coord.resize(2);
        double x = std::get<1>(coordinates_gf[0])(i);
        double y = std::get<1>(coordinates_gf[1])(i);
        double z = std::get<1>(coordinates_gf[2])(i);
        double r = std::sqrt(x * x + y * y);
        coord[0] = r;
        coord[1] = z;
      } else if constexpr (CS == CoordinateSystem::Radial) {
        coord.resize(1);
        double x = std::get<1>(coordinates_gf[0])(i);
        double y = std::get<1>(coordinates_gf[1])(i);
        double r = std::sqrt(x * x + y * y);
        coord[0] = r;
      }
      auto fun = funcs[u];

      unk(i) = fun(coord, time);
    }
  }
}

/**
 * @brief
 *
 * @tparam CS
 */
template <CoordinateSystem CS>
void FittedLoading<CS>::initialize() {}
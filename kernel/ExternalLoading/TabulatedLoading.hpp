/**
 * @file TabulatedLoading.hpp
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
#include "Inputs/HDF54Sloth.hpp"
#include "Interpolators/MultiLinearInterpolator.hpp"
#include "Parameters/Parameter.hpp"

#pragma once

template <CoordinateSystem CS>
class TabulatedLoading {
 public:
  std::vector<std::function<double(const std::vector<double>&, double)>> funcs;
  std::vector<FlattenedTensor<double>> vector_tensors;
  // liste des variables un pour chaque, si pas utile pas rempli
  std::vector<std::vector<double>> vector_time;
  std::vector<std::vector<double>> vector_xr;
  std::vector<std::vector<double>> vector_y;
  std::vector<std::vector<double>> vector_z;

  std::vector<std::string> hdf5_filename;
  std::vector<std::string> hdf5_datasetname;
  const Parameters& params_;
  void do_time_step(double& time, std::vector<mfem::Vector>& unknown_,
                    const std::vector<std::vector<std::string>>& unks_info,
                    const std::vector<std::tuple<std::string, mfem::Vector>>& coordinates_gf);
  void initialize();
  void get_parameters();

  TabulatedLoading(const Parameters& params);
};

/**
 * @brief Construct a new Fitted Loading<CS>:: Fitted Loading object
 *
 * @tparam CS
 * @param params
 */
template <CoordinateSystem CS>
TabulatedLoading<CS>::TabulatedLoading(const Parameters& params) : params_(params) {
  this->get_parameters();
}

/**
 * @brief
 *
 * @tparam CS
 */
template <CoordinateSystem CS>
void TabulatedLoading<CS>::get_parameters() {
  this->funcs = this->params_.template get_param_value_or_default<
      std::vector<std::function<double(const std::vector<double>&, double)>>>("LoadingFunctions",
                                                                              {});
  this->hdf5_filename = this->params_.template get_param_value_or_default<std::vector<std::string>>(
      "hdf5_filename", {});
  this->hdf5_datasetname =
      this->params_.template get_param_value_or_default<std::vector<std::string>>(
          "hdf5_datasetname", {});
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
void TabulatedLoading<CS>::do_time_step(
    double& time, std::vector<mfem::Vector>& unknown_,
    const std::vector<std::vector<std::string>>& unks_info,
    const std::vector<std::tuple<std::string, mfem::Vector>>& coordinates_gf) {
  std::vector<double> coord;
  std::vector<std::size_t> lower_indices;
  std::vector<std::vector<double>> grid_values;
  std::vector<double> point_to_interpolate;
  for (std::size_t u = 0; u < unknown_.size(); u++) {  // loop on unknown
    auto& unk = unknown_[u];
    for (std::size_t i = 0; i < unk.Size(); i++) {  // loop on each point of the unknows

      std::size_t dim = coordinates_gf.size();  // problem dimension

      if constexpr (CS == CoordinateSystem::Cartesian) {
        // coord.resize(dim);
        // for (std::size_t j = 0; j < dim; j++) {
        //   coord[j] = std::get<1>(coordinates_gf[j])(i);  // get physical coordinates (x,y,z)
        // }
        throw std::runtime_error("CoordinateSystem::Cartesian not available");
      } else if constexpr (CS == CoordinateSystem::Axisymmetric) {
        // coord.resize(2);
        // double x = std::get<1>(coordinates_gf[0])(i);
        // double y = std::get<1>(coordinates_gf[1])(i);
        // double z = std::get<1>(coordinates_gf[2])(i);
        // double r = std::sqrt(x * x + y * y);
        // coord[0] = r;
        // coord[1] = z;
        throw std::runtime_error("CoordinateSystem::Cartesian not available");
      } else if constexpr (CS == CoordinateSystem::Radial) {
        coord.resize(1);
        lower_indices.resize(2);
        grid_values.resize(2);
        point_to_interpolate.resize(2);
        double x = std::get<1>(coordinates_gf[0])(i);
        double y = std::get<1>(coordinates_gf[1])(i);
        double r = std::sqrt(x * x + y * y);
        coord[0] = r;
        auto ftensor = vector_tensors[u];
        auto vect_time = vector_time[u];
        auto vect_x = vector_xr[u];

        auto it_t = std::lower_bound(vect_time.begin(), vect_time.end(), time);
        std::size_t index_t = it_t - vect_time.begin();
        index_t = std::max(std::size_t(0),
                           std::min(index_t == 0 ? 0 : index_t - 1, vect_time.size() - 2));

        auto it_r = std::lower_bound(vect_x.begin(), vect_x.end(), r);
        std::size_t index_r = it_r - vect_x.begin();
        index_r =
            std::max(std::size_t(0), std::min(index_r == 0 ? 0 : index_r - 1, vect_x.size() - 2));

        lower_indices[1] = index_t;
        lower_indices[0] = index_r;
        grid_values[1] = vect_time;
        grid_values[0] = vect_x;
        point_to_interpolate[1] = time;
        point_to_interpolate[0] = r;
        std::vector<double> alpha =
            MultiLinearInterpolator<std::monostate>::computeInterpolationCoefficients(
                2, point_to_interpolate, lower_indices, grid_values);

        // unk(i) = unk(i) * std::pow(MultiLinearInterpolator<std::monostate>::computeInterpolation(
        //                       2, lower_indices, alpha, ftensor),0.1);
              unk(i) = MultiLinearInterpolator<std::monostate>::computeInterpolation(
                              2, lower_indices, alpha, ftensor);
                              // std::cout << unk(i) << std::endl;
      }

      // unk(i) = 0.;
    }
  }
}

/**
 * @brief
 *
 * @tparam CS
 */
template <CoordinateSystem CS>
void TabulatedLoading<CS>::initialize() {
  std::size_t var_num = this->hdf5_filename.size();
  vector_tensors.resize(var_num);
  vector_time.resize(var_num);
  vector_xr.resize(var_num);
  vector_y.resize(var_num);
  vector_z.resize(var_num);
  for (std::size_t i = 0; i < this->hdf5_filename.size(); i++) {
    HDF54Sloth<std::monostate>::get_data_from_HDF5(this->hdf5_filename[i],
                                                   this->hdf5_datasetname[i], vector_tensors[i]);
    HDF54Sloth<std::monostate>::get_data_from_HDF5(this->hdf5_filename[i], "t", vector_time[i]);
    if constexpr (CS == CoordinateSystem::Cartesian) {
      throw std::runtime_error("CoordinateSystem::Cartesian not available");
    } else if constexpr (CS == CoordinateSystem::Axisymmetric) {
      throw std::runtime_error("CoordinateSystem::Cartesian not available");
      // HDF54Sloth<std::monostate>::get_data_from_HDF5(this->hdf5_filename[i], "x", vector_xr[i]);
      // HDF54Sloth<std::monostate>::get_data_from_HDF5(this->hdf5_filename[i], "y", vector_y[i]);
      // HDF54Sloth<std::monostate>::get_data_from_HDF5(this->hdf5_filename[i], "z", vector_z[i]);
    } else if constexpr (CS == CoordinateSystem::Radial) {
      HDF54Sloth<std::monostate>::get_data_from_HDF5(this->hdf5_filename[i], "r", vector_xr[i]);
    }
  }
}
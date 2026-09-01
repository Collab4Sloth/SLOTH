/**
 * @file postprocessing.tpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief  Post-processing features
 * @version 0.1
 * @date 2025-09-05
 *
 * @copyright CEA (C) 2025
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

#pragma once
#include <filesystem>  // NOLINT [avoid  <filesystem> is an unapproved C++17 header.]
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "Spatial/Spatial.hpp"
#include "Utils/Utils.hpp"
#include "Variables/Variable.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new Post Processing<T,DC,DIM>:: Post Processing object
 *
 * @tparam T mfem FECollection
 * @tparam DC mfem DataCollection
 * @tparam DIM Spatial dimension
 * @param space Collection of SpatialDiscretization objects associated with the Variables of the
 * Problem
 * @param params Paramters used by the PostProcessing object
 */
template <class T, class DC, int DIM>
PostProcessing<T, DC, DIM>::PostProcessing(SpatialDiscretization<T, DIM>* space,
                                           const Parameters& params)
    : params_(params) {
  this->dc_ = std::make_shared<DC>(params.get_param_value<std::string>("calculation_path"),
                                   space->get_mesh());
  this->get_parameters();

  this->clean_output_directory();

  this->dc_->SetPrefixPath(this->main_folder_path_);
  this->dc_->SetLevelsOfDetail(this->level_of_detail_);
  this->dc_->SetDataFormat(mfem::VTKFormat::BINARY);
  this->dc_->SetHighOrderOutput(true);
  this->post_processing_directory_ = this->main_folder_path_ + "/" + this->calculation_path_;

  if (mfem::Mpi::WorldRank() == 0) {
    std::filesystem::path dir_path = std::filesystem::path(this->post_processing_directory_);
    if (!std::filesystem::exists(dir_path)) {
      std::filesystem::create_directories(dir_path);
    }
  }
  // Wait creation of the directory before continue
  MPI_Barrier(MPI_COMM_WORLD);
}

/**
 * @brief Get specialized parameters
 *
 * @tparam T mfem FECollection
 * @tparam DC mfem DataCollection
 * @tparam DIM Spatial dimension
 */
template <class T, class DC, int DIM>
void PostProcessing<T, DC, DIM>::get_specialized_parameters() {
  this->enable_compute_energies_ =
      this->params_.template get_param_value_or_default<bool>("enable_compute_energies", true);

  if (this->params_.has_parameter("iso_val_to_compute")) {
    this->iso_val_to_compute_ =
        this->params_.template get_param_value<MapStringDouble>("iso_val_to_compute");
  }

  if (this->params_.has_parameter("integral_to_compute")) {
    this->integral_to_compute_ =
        this->params_.template get_param_value<MapString2Double>("integral_to_compute");
    for (const auto& [variable, bounds] : this->integral_to_compute_) {
      const auto& [lower_bound, upper_bound] = bounds;
      std::string error_msg =
          "Error with variable " + variable + ": Lower bound must lower than Upper bound";
      MFEM_VERIFY(upper_bound >= lower_bound, error_msg.c_str());
    }
  }
}
/**
 * @brief Get restart parameters
 *
 * @tparam T mfem FECollection
 * @tparam DC mfem DataCollection
 * @tparam DIM Spatial dimension
 */
template <class T, class DC, int DIM>
void PostProcessing<T, DC, DIM>::get_restart_parameters() {
  if (this->params_.has_parameter("iterations_list_save_gf")) {
    this->iterations_list_save_gf_ =
        this->params_.template get_param_value<vInt>("iterations_list_save_gf");

    this->gf_folder_path_ =
        this->params_.template get_param_value_or_default<std::string>("gf_folder_path", "GF");

    std::filesystem::path dir = this->gf_folder_path_;
    if (!std::filesystem::exists(dir)) {
      std::filesystem::create_directories(dir);
    }
  }
}

/**
 * @brief Get common parameters
 *
 * @tparam T mfem FECollection
 * @tparam DC mfem DataCollection
 * @tparam DIM Spatial dimension
 */
template <class T, class DC, int DIM>
void PostProcessing<T, DC, DIM>::get_common_parameters() {
  this->main_folder_path_ = this->params_.template get_param_value<std::string>("main_folder_path");
  this->calculation_path_ = this->params_.template get_param_value<std::string>("calculation_path");

  this->level_of_detail_ =
      this->params_.template get_param_value_or_default<int>("level_of_detail", 1);
  this->force_clean_output_dir_ =
      this->params_.template get_param_value_or_default<bool>("force_clean_output_dir", false);

  this->enable_save_specialized_at_iter_ = this->params_.template get_param_value_or_default<bool>(
      "enable_save_specialized_at_iter", true);
  std::string error_msg =
      "At least one the following parameter is expected: frequency, iterations_list, times_list";
  MFEM_VERIFY(this->params_.has_parameter("frequency") ||
                  this->params_.has_parameter("iterations_list") ||
                  this->params_.has_parameter("times_list"),
              error_msg.c_str());
  if (this->params_.has_parameter("frequency")) {
    this->frequency_ = this->params_.template get_param_value<int>("frequency");
  } else {
    this->frequency_ = 1;
    if (this->params_.has_parameter("iterations_list")) {
      this->iterations_list_ = this->params_.template get_param_value<vInt>("iterations_list");
    }
    if (this->params_.has_parameter("times_list")) {
      this->times_list_ = this->params_.template get_param_value<vDouble>("times_list");
    }
  }
}

/**
 * @brief Get parameters for PostProcessing object
 *
 * @tparam T mfem FECollection
 * @tparam DC mfem DataCollection
 * @tparam DIM Spatial dimension
 */
template <class T, class DC, int DIM>
void PostProcessing<T, DC, DIM>::get_parameters() {
  // Global
  this->get_common_parameters();

  // Specialized (CSV)
  this->get_specialized_parameters();

  // Restart
  this->get_restart_parameters();
}

/**
 * @brief Save Variables at a given iteration/time
 *
 * @tparam T mfem FECollection
 * @tparam DC mfem DataCollection
 * @tparam DIM Spatial dimension
 * @param vars
 * @param iter
 * @param time
 */
template <class T, class DC, int DIM>
void PostProcessing<T, DC, DIM>::save_variables(Variables<T, DIM>& vars, const int& iter,
                                                const double& time) {
  // VTK
  if (this->need_to_be_saved(iter, time)) {
    this->dc_->SetCycle(iter);
    this->dc_->SetTime(time);
    std::map<std::string, mfem::ParGridFunction*> map_var = vars.get_map_gridfunction();
    for (auto& [name, gf_ptr] : map_var) {
      this->dc_->RegisterField(name, gf_ptr);
    }
    this->dc_->Save();
  }
  // Restart
  if (std::ranges::find(this->iterations_list_save_gf_, iter) !=
      this->iterations_list_save_gf_.end()) {
    std::map<std::string, mfem::ParGridFunction*> map_var = vars.get_map_gridfunction();
    std::filesystem::path output_dir_path(this->gf_folder_path_);
    for (auto& [name, gf_ptr] : map_var) {
      std::filesystem::path filename = name + "_" + std::to_string(iter) + ".gf";
      std::filesystem::path full_path = output_dir_path / filename;
      gf_ptr->Save(full_path.string().c_str(), 16);
    }
  }
}

/**
 * @brief Get the frequency of post-processing in terms of number of iterations (1 means each
 * iteration)
 *
 * @tparam T mfem FECollection
 * @tparam DC mfem DataCollection
 * @tparam DIM Spatial dimension
 * @return int The frequency of post-processing
 */
template <class T, class DC, int DIM>
int PostProcessing<T, DC, DIM>::get_frequency() {
  return this->frequency_;
}

/**
 * @brief Get the isovalues to compute
 *
 * @tparam T mfem FECollection
 * @tparam DC mfem DataCollection
 * @tparam DIM Spatial dimension
 * @return std::map<std::string, double>
 */
template <class T, class DC, int DIM>
std::map<std::string, double> PostProcessing<T, DC, DIM>::get_iso_val_to_compute() {
  return this->iso_val_to_compute_;
}

/**
 * @brief Get the integrals to compute over the domain
 *
 * @tparam T mfem FECollection
 * @tparam DC mfem DataCollection
 * @tparam DIM Spatial dimension
 * @return std::map<std::string, double>
 */
template <class T, class DC, int DIM>
std::map<std::string, std::tuple<double, double>>
PostProcessing<T, DC, DIM>::get_integral_to_compute() {
  return this->integral_to_compute_;
}

/**
 * @brief Return the post-processing directory
 *
 * @tparam T mfem FECollection
 * @tparam DC mfem DataCollection
 * @tparam DIM Spatial dimension
 * @return std::string
 */
template <class T, class DC, int DIM>
std::string PostProcessing<T, DC, DIM>::get_post_processing_directory() {
  return this->post_processing_directory_;
}

/**
 * @brief Indicate if specialized values must be saved at each iteration or not
 *
 * @tparam T mfem FECollection
 * @tparam DC mfem DataCollection
 * @tparam DIM Spatial dimension
 * @return true
 * @return false
 */
template <class T, class DC, int DIM>
bool PostProcessing<T, DC, DIM>::get_enable_save_specialized_at_iter() {
  return this->enable_save_specialized_at_iter_;
}

/**
 * @brief Indicate if energies must be calculated
 *
 * @tparam T mfem FECollection
 * @tparam DC mfem DataCollection
 * @tparam DIM Spatial dimension
 * @return true
 * @return false
 */
template <class T, class DC, int DIM>
bool PostProcessing<T, DC, DIM>::get_enable_compute_energies() {
  return this->enable_compute_energies_;
}

/**
 * @brief Check if results have to be saved at iteration
 *
 * @tparam T mfem FECollection
 * @tparam DC mfem DataCollection
 * @tparam DIM Spatial dimension
 * @param iteration The current iteration.
 * @return true
 * @return false
 */
template <class T, class DC, int DIM>
bool PostProcessing<T, DC, DIM>::need_to_be_saved(const int iteration, const double time) {
  // Tolerance used to detect time to save
  const double epsilon = 1.e-12;

  bool check_frequency = (iteration % this->frequency_ == 0);
  bool check_iterations_list =
      std::ranges::find(this->iterations_list_, iteration) != this->iterations_list_.end();
  bool check_times_list =
      std::any_of(this->times_list_.begin(), this->times_list_.end(),
                  [time, epsilon](double t) { return std::abs(t - time) < epsilon; });
  bool check = check_frequency;
  if (!this->iterations_list_.empty() || !this->times_list_.empty()) {
    check &= (check_iterations_list || check_times_list);
  }
  return check;
}

/**
 * @brief Export specialized results in CSV files
 *
 * @tparam T mfem FECollection
 * @tparam DC mfem DataCollection
 * @tparam DIM Spatial dimension
 * @param mmap_results
 * @param filename
 */
template <class T, class DC, int DIM>
void PostProcessing<T, DC, DIM>::save_specialized(
    const std::multimap<IterationKey, SpecializedValue>& mmap_results, std::string filename) {
  std::filesystem::path file = std::filesystem::path(this->post_processing_directory_) / filename;

  if (!std::filesystem::exists(file)) {
    std::ostringstream text2fic;
    // File doesn't exist
    std::ofstream fic(file, std::ios::out);

    if (fic.is_open()) {
      ////////////////////////////////////////////
      // Headers
      ////////////////////////////////////////////
      auto key0 = mmap_results.begin()->first;
      auto value0 = mmap_results.begin()->second;
      auto range0 = mmap_results.equal_range(key0);
      text2fic << key0.iter_.first << "," << key0.time_step_.first << "," << key0.time_.first;
      for (auto it = range0.first; it != range0.second; ++it) {
        const auto& value = it->second;
        text2fic << "," << value.first;
      }
      text2fic << "\n";
      fic << text2fic.str();
      fic.close();
    }
  }
  // File already exists
  std::ofstream fic(file, std::ios::out | std::ios::app);

  if (!fic.is_open()) {
    std::string msg = "Unable to open file: " + filename;
    mfem::mfem_error(msg.c_str());
  }
  ////////////////////////////////////////////
  // Values
  ////////////////////////////////////////////
  std::ostringstream text2fic;
  std::set<IterationKey> already_seen_keys;
  for (const auto& [key, value] : mmap_results) {
    if (already_seen_keys.find(key) != already_seen_keys.end()) {
      continue;
    }
    auto range = mmap_results.equal_range(key);
    text2fic << key.iter_.second << "," << key.time_step_.second << "," << key.time_.second;
    for (auto it = range.first; it != range.second; ++it) {
      const auto& value = it->second;
      text2fic << "," << value.second;
    }
    text2fic << "\n";
    already_seen_keys.insert(key);
  }
  fic << text2fic.str();
  fic.close();
}

/**
 * @brief Clean output_directory before calculation
 *
 * @tparam T mfem FECollection
 * @tparam DC mfem DataCollection
 * @tparam DIM Spatial dimension
 */
template <class T, class DC, int DIM>
void PostProcessing<T, DC, DIM>::clean_output_directory() {
  int rank = mfem::Mpi::WorldRank();
  if (rank == 0) {
    if (this->force_clean_output_dir_) {
      auto output_dir_path = std::filesystem::path(this->main_folder_path_);
      std::error_code ec;
      std::filesystem::remove_all(output_dir_path, ec);
      if (ec) {
        auto msg = ec.message();
        mfem::mfem_error(msg.c_str());
      }
    }
  }
}

/**
 * @brief Collect this Problem's grid functions into a shared field map.
 *
 * @tparam T mfem FECollection
 * @tparam DC mfem DataCollection
 * @tparam DIM Spatial dimension
 * @param vars Variables to collect fields from.
 * @param all_fields Output map, accumulated across multiple Problems for
 *                   a unified VTK save (not cleared by this call).
 */
template <class T, class DC, int DIM>
void PostProcessing<T, DC, DIM>::collect_vtk_fields(
    Variables<T, DIM>& vars, std::map<std::string, mfem::ParGridFunction*>& all_fields) {
  auto map_var = vars.get_map_gridfunction();
  all_fields.insert(map_var.begin(), map_var.end());
}

/**
 * @brief Return the shared DataCollection used for VTK output.
 *
 * @tparam T mfem FECollection
 * @tparam DC mfem DataCollection
 * @tparam DIM Spatial dimension
 * @return std::shared_ptr<DC> The DataCollection instance.
 */
template <class T, class DC, int DIM>
std::shared_ptr<DC> PostProcessing<T, DC, DIM>::get_shared_dc() {
  return this->dc_;
}
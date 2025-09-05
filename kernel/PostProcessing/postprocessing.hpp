/**
 * @file postprocessing.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief  Post-processing features
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
 * @brief Class used to manage post-processing of SLOTH simulations
 *
 * @tparam T mfem FECollection
 * @tparam DC mfem DataCollection
 * @tparam DIM Spatial dimension
 */
template <class T, class DC, int DIM>
class PostProcessing : public DC {
 private:
  std::string main_folder_path_;
  std::string calculation_path_;

  int frequency_;
  std::vector<int> iterations_list_;
  std::vector<double> times_list_;

  int level_of_detail_;
  bool enable_save_specialized_at_iter_;
  bool force_clean_output_dir_;
  std::map<std::string, double> iso_val_to_compute_;
  std::map<std::string, std::tuple<double, double>> integral_to_compute_;

  const Parameters& params_;
  std::map<std::string, mfem::ParGridFunction> fields_to_save_;
  std::string post_processing_directory_;
  void get_parameters();

  bool need_to_be_saved(const int iteration, const double time);

  void clean_output_directory();

 public:
  // Explicitly default the move constructor and move assignment operator
  // Usefull to define several PostProcessing objects in std::vector over a loop
  PostProcessing(PostProcessing&&) = default;
  PostProcessing& operator=(PostProcessing&&) = default;

  PostProcessing(SpatialDiscretization<T, DIM>* space, const Parameters& params);
  void save_variables(const Variables<T, DIM>& vars, const int& iter, const double& time);
  void save_specialized(const std::multimap<IterationKey, SpecializedValue>& mmap_results,
                        std::string filename = "time_specialized.csv");
  void save_iso_specialized(const std::multimap<IterationKey, SpecializedValue>& mmap_results,
                            std::string filename = "iso.csv");
  int get_frequency();
  std::string get_post_processing_directory();
  bool get_enable_save_specialized_at_iter();
  std::map<std::string, double> get_iso_val_to_compute();
  std::map<std::string, std::tuple<double, double>> get_integral_to_compute();

  ~PostProcessing();
};

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

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
    : DC(params.get_param_value<std::string>("calculation_path"), space->get_mesh()),
      params_(params) {
  this->get_parameters();

  this->clean_output_directory();

  this->SetPrefixPath(this->main_folder_path_);
  this->SetLevelsOfDetail(this->level_of_detail_);
  this->SetDataFormat(mfem::VTKFormat::BINARY);
  this->SetHighOrderOutput(true);
  this->post_processing_directory_ = this->main_folder_path_ + "/" + this->calculation_path_;
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
  this->main_folder_path_ = this->params_.template get_param_value<std::string>("main_folder_path");
  this->calculation_path_ = this->params_.template get_param_value<std::string>("calculation_path");

  this->level_of_detail_ =
      this->params_.template get_param_value_or_default<int>("level_of_detail", 1);
  this->enable_save_specialized_at_iter_ = this->params_.template get_param_value_or_default<bool>(
      "enable_save_specialized_at_iter", false);
  this->force_clean_output_dir_ =
      this->params_.template get_param_value_or_default<bool>("force_clean_output_dir", false);
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
void PostProcessing<T, DC, DIM>::save_variables(const Variables<T, DIM>& vars, const int& iter,
                                                const double& time) {
  if (this->need_to_be_saved(iter, time)) {
    this->SetCycle(iter);
    this->SetTime(time);
    std::map<std::string, mfem::ParGridFunction> map_var = vars.get_map_gridfunction();
    for (auto& [name, gf] : map_var) {
      this->RegisterField(name, &gf);
    }
    this->Save();
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
 * @brief Destroy the Post Processing:: Post Processing object
 *
 */
template <class T, class DC, int DIM>
PostProcessing<T, DC, DIM>::~PostProcessing() {}

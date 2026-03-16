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

#include "Parameters/Parameters.hpp"
#include "Spatial/Spatial.hpp"
#include "Utils/Utils.hpp"
#include "Variables/Variable.hpp"
#include "Variables/Variables.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Class used to manage post-processing of SLOTH simulations
 *
 * @tparam T mfem FECollection
 * @tparam DC mfem DataCollection
 * @tparam DIM Spatial dimension
 */
template <class T, class DC, int DIM>
class PostProcessing {
 private:
  std::shared_ptr<DC> dc_;
  std::string main_folder_path_;
  std::string calculation_path_;

  int frequency_;
  std::vector<int> iterations_list_;
  std::vector<double> times_list_;

  int level_of_detail_;
  bool enable_compute_energies_;
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
  bool get_enable_compute_energies();
  std::map<std::string, double> get_iso_val_to_compute();
  std::map<std::string, std::tuple<double, double>> get_integral_to_compute();

  virtual ~PostProcessing() = default;
};

#include "PostProcessing/postprocessing.tpp"

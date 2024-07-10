/*
 * Copyright © CEA 2023
 *
 * \brief Post-processing used by phase-field models
 *
 * \file PostProcessing.hpp
 * \author ci230846
 * \date 27/03/2023
 */

#pragma once
#include <filesystem>
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

#include "Spatial/Spatial.hpp"
#include "Utils/PhaseFieldOptions.hpp"
#include "Variables/Variable.hpp"
#include "mfem.hpp" // NOLINT [no include the directory when naming mfem include file]

template <class T, class DC, int DIM>
class PostProcessing : public DC {
 private:
  std::map<std::string, mfem::GridFunction> fields_to_save_;
  int frequency_;
  std::string post_processing_directory_;
  bool need_to_be_saved(const int& iteration);

 public:
  PostProcessing(const std::string& main_folder_path, const std::string& calculation_path,
                 SpatialDiscretization<T, DIM>* space, const int& frequency,
                 const int& level_of_detail);
  void save_variables(const Variables<T, DIM>& vars, const int& iter, const double& time);
  void save_specialized(const std::multimap<IterationKey, SpecializedValue>& mmap_results);
  int get_frequency();

  ~PostProcessing();
};

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

/**
 * @brief Construct a new Post Processing:: Post Processing object
 *
 * @param main_folder_path
 * @param calculation_path
 * @param mesh
 * @param level_of_detail
 */
template <class T, class DC, int DIM>
PostProcessing<T, DC, DIM>::PostProcessing(const std::string& main_folder_path,
                                           const std::string& calculation_path,
                                           SpatialDiscretization<T, DIM>* space,
                                           const int& frequency, const int& level_of_detail)
    : DC(calculation_path, &space->get_mesh()), frequency_(frequency) {
  this->SetPrefixPath(main_folder_path);
  this->SetLevelsOfDetail(level_of_detail);
  this->SetDataFormat(mfem::VTKFormat::BINARY);
  this->SetHighOrderOutput(true);
  this->post_processing_directory_ = main_folder_path + "/" + calculation_path;
}

/**
 * @brief save variables objet at given iter/time
 *
 * @param vars
 * @param iter
 * @param time
 */
template <class T, class DC, int DIM>
void PostProcessing<T, DC, DIM>::save_variables(const Variables<T, DIM>& vars, const int& iter,
                                                const double& time) {
  if (this->need_to_be_saved(iter)) {
    this->SetCycle(iter);
    this->SetTime(time);
    auto map_var = vars.get_map_gridfunction();
    for (auto [name, gf] : map_var) {
      this->RegisterField(name, &gf);
      this->Save();
    }
  }
}

/**
 * @brief Get the frequency of post-processing in terms of number of iterations (1 means each
 * iteration)
 *
 * @return int
 */
template <class T, class DC, int DIM>
int PostProcessing<T, DC, DIM>::get_frequency() {
  return this->frequency_;
}

/**
 * @brief check if results have to be saved at iteration
 *
 * @param iteration
 * @return true
 * @return false
 */
template <class T, class DC, int DIM>
bool PostProcessing<T, DC, DIM>::need_to_be_saved(const int& iteration) {
  bool check = (iteration % this->frequency_ == 0);
  return check;
}

/**
 * @brief Export specialized results in CSV files
 *
 * @tparam T
 * @tparam DC
 * @tparam DIM
 * @param filename
 * @param tup
 */
template <class T, class DC, int DIM>
void PostProcessing<T, DC, DIM>::save_specialized(
    const std::multimap<IterationKey, SpecializedValue>& mmap_results) {
  const std::string& filename = "time_specialized.csv";
  std::filesystem::path file = std::filesystem::path(this->post_processing_directory_) / filename;
  std::ofstream fic(file, std::ios::out | std::ios::trunc);

  if (!fic.is_open()) {
    throw std::runtime_error("Unable to open file: " + filename);
  }

  std::ostringstream text2fic;
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
  ////////////////////////////////////////////
  // Values
  ////////////////////////////////////////////
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
}

/**
 * @brief Destroy the Post Processing:: Post Processing object
 *
 */
template <class T, class DC, int DIM>
PostProcessing<T, DC, DIM>::~PostProcessing() {}

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
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include "../Spatial/Spatial.hpp"
#include "../Utils/PhaseFieldOptions.hpp"
#include "../Variables/Variable.hpp"
#include "../Variables/Variables.hpp"
#include "mfem.hpp"

template <class T, class DC, int DIM>
class PostProcessing : public DC {
 private:
  std::map<std::string, mfem::GridFunction> fields_to_save_;
  int frequency_;
  std::string post_processing_directory_;

 public:
  PostProcessing(const std::string& main_folder_path, const std::string& calculation_path,
                 SpatialDiscretization<T, DIM>* space, const int& frequency,
                 const int& level_of_detail);
  void save_variables(const Variables<T, DIM>& vars, const int& iter, const double& time);
  void create_csv(const std::string& filename, const std::string& headers);
  void export_csv(const std::string& filename,
                  const std::map<std::tuple<int, double, double>, double>& map_results);
  int get_frequency();
  bool need_to_be_saved(const int& iteration);
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
  this->SetCycle(iter);
  this->SetTime(time);
  auto map_var = vars.get_map_gridfunction();
  for (auto [name, gf] : map_var) {
    this->RegisterField(name, &gf);
    this->Save();
  }
}

/**
 * @brief create csv file from name, overwrite if already exist
 *
 * @tparam T
 * @tparam DC
 * @tparam DIM
 * @param filename
 */
template <class T, class DC, int DIM>
void PostProcessing<T, DC, DIM>::create_csv(const std::string& filename,
                                            const std::string& headers) {
  std::ios_base::openmode mode = std::ios::out | std::ios::trunc;
  std::ofstream fic;
  fic.open(this->post_processing_directory_ + "/" + filename, mode);
  fic << headers;
  fic << std::endl;
  fic.close();
}

/**
 * @brief export results into CSV file
 *
 * @tparam T
 * @tparam DC
 * @tparam DIM
 * @param map_results
 */
template <class T, class DC, int DIM>
void PostProcessing<T, DC, DIM>::export_csv(
    const std::string& filename,
    const std::map<std::tuple<int, double, double>, double>& map_results) {
  std::ios_base::openmode mode = std::ios::out | std::ios::app;
  std::ofstream fic;
  fic.open(this->post_processing_directory_ + "/" + filename, mode);
  std::ostringstream text2fic;
  // TODO(CCI) : template + forward ?
  for (auto [key, value] : map_results) {
    text2fic << std::get<0>(key) << " " << std::get<1>(key) << " " << std::get<2>(key) << " "
             << value << "\n";
  }
  fic << text2fic.str();
  fic.close();
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
 * @brief Destroy the Post Processing:: Post Processing object
 *
 */
template <class T, class DC, int DIM>
PostProcessing<T, DC, DIM>::~PostProcessing() {}

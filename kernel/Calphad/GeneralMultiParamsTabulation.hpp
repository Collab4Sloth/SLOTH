/**
 * @file GeneralMultiParamsTabulation.hpp
 * @author cp273896  (clement.plumecocq@cea.fr)
 * @brief Tabulated thermodynamic description for an binary isothermal system
 * @version 0.1
 * @date 2025-03-17
 *
 * Copyright CEA (c) 2025
 *
 */
#include <H5Cpp.h>

#include <algorithm>
#include <boost/multi_array.hpp>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Calphad/CalphadBase.hpp"
#include "Calphad/CalphadUtils.hpp"
#include "Inputs/HDF54Sloth.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"

#pragma once

template <typename T>
class GeneralMultiParamsTabulation : CalphadBase<T> {
 private:
  std::vector<double> tab_T, tab_x1, tab_x2;

  boost::multi_array<double, 3> array_g;
  std::map<std::string, boost::multi_array<double, 3>> array_mu;
  std::map<std::string, boost::multi_array<double, 3>> array_mobility;

  std::unique_ptr<CalphadUtils<T>> CU_;
  void compute(
      const size_t nb_nodes, const std::vector<T>& tp_gf,
      const std::vector<std::tuple<std::string, std::string>>& chemical_system,
      std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system);

  void check_variables_consistency(
      std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system);

 public:
  explicit GeneralMultiParamsTabulation(const Parameters& params);

  void initialize() override;

  void execute(const int dt, const std::vector<T>& aux_gf,
               const std::vector<std::tuple<std::string, std::string>>& chemical_system,
               std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>&
                   output_system) override;

  void finalize() override;

  ////////////////////////////////

  void get_parameters() override;

  ////////////////////////////////

  ~GeneralMultiParamsTabulation();
};

////////////////////////////////
////////////////////////////////

/**
 * @brief Get the parameters associated with the GeneralMultiParamsTabulation object
 *
 * @tparam T
 */
template <typename T>
void GeneralMultiParamsTabulation<T>::get_parameters() {
  this->description_ = this->params_.template get_param_value_or_default<std::string>(
      "description",
      "Analytical thermodynamic description for an ideal solution using tabulated data ");
}

////////////////////////////////
////////////////////////////////
/**
 * @brief Construct a new GeneralMultiParamsTabulation::GeneralMultiParamsTabulation object
 *
 * @param params
 */
template <typename T>
GeneralMultiParamsTabulation<T>::GeneralMultiParamsTabulation(const Parameters& params) : CalphadBase<T>(params) {
  this->CU_ = std::make_unique<CalphadUtils<T>>();
  this->get_parameters();
}

/**
 * @brief Initialization of the thermodynamic calculation
 *
 * @tparam T
 */
template <typename T>
void GeneralMultiParamsTabulation<T>::initialize() {
  Catch_Time_Section("GeneralMultiParamsTabulation::initialize");
  std::string filename = this->params_.template get_param_value_or_default<std::string>(
      "data filename", "no input file");

  H5std_string DATASET_NAME;
  const H5std_string FILE_NAME(filename);

  HDF54Sloth<3> hdf5_for_2D_table;

  std::vector<std::string> list_of_element = {"O", "PU", "U"};
  std::vector<std::string> list_of_dataset_name_mu = {"mu_O", "mu_PU", "mu_U"};
  std::vector<std::string> list_of_dataset_name_mob = {"Mo", "Mpu", "Mu"}; 

  for (size_t i = 0; i < list_of_element.size(); i++) {
    boost::multi_array<double, 3> myArray;

    DATASET_NAME = list_of_dataset_name_mu[i];
    hdf5_for_2D_table.get_data_from_HDF5(FILE_NAME, DATASET_NAME, myArray);
    array_mu.emplace(list_of_element[i], myArray);

    DATASET_NAME = list_of_dataset_name_mob[i];
    hdf5_for_2D_table.get_data_from_HDF5(FILE_NAME, DATASET_NAME, myArray);
    array_mobility.emplace(list_of_element[i], myArray);
  }

  boost::multi_array<double, 3> array_muO;
  DATASET_NAME = "G";
  hdf5_for_2D_table.get_data_from_HDF5(FILE_NAME, DATASET_NAME, array_g);

  HDF54Sloth<1> hdf5_for_vector;
  DATASET_NAME = "xO";
  hdf5_for_vector.get_data_from_HDF5(FILE_NAME, DATASET_NAME, tab_x1);
  DATASET_NAME = "xU";
  hdf5_for_vector.get_data_from_HDF5(FILE_NAME, DATASET_NAME, tab_x2);
  DATASET_NAME = "T";
  hdf5_for_vector.get_data_from_HDF5(FILE_NAME, DATASET_NAME, tab_T);
}

/**
 * @brief Main method to calculate equilibrium states
 *
 * @tparam T
 * @param dt
 * @param aux_gf
 * @param chemical_system
 * @param output_system
 */
template <typename T>
void GeneralMultiParamsTabulation<T>::execute(
    const int dt, const std::vector<T>& tp_gf,
    const std::vector<std::tuple<std::string, std::string>>& chemical_system,
    std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system) {
  // Clear containers and recalculation of the numbers of nodes
  const size_t nb_nodes = this->CU_->get_size(tp_gf[0]);
  this->clear_containers();
  if (dt == 1) {
    this->check_variables_consistency(output_system);
  }
  // Thermodynamic Calculations
  this->compute(nb_nodes, tp_gf, chemical_system, output_system);

  // Use containers to update output_system
  this->update_outputs(nb_nodes, output_system);
}

/**
 * @brief Compute the CALPHAD contributions
 *
 * @tparam T
 * @param tp_gf
 * @param chemical_system
 * @param output_system
 */
template <typename T>
void GeneralMultiParamsTabulation<T>::compute(
    const size_t nb_nodes, const std::vector<T>& tp_gf,
    const std::vector<std::tuple<std::string, std::string>>& chemical_system,
    std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system) {
  Catch_Time_Section("GeneralMultiParamsTabulation::compute");

  // Let us assume an ideal mixing solution
  const std::string& phase = "C1_MO2";
  const std::vector<std::string> energy_names = {"G"};
  std::vector<double> tp_gf_at_node(tp_gf.size());

  std::vector<int> sorted_n_t_p = this->CU_->sort_nodes(tp_gf[0], tp_gf[1], "No", "No");

  // Process CALPHAD calculations for each node
  for (const auto& id : sorted_n_t_p) {
    // Populate tp_gf_at_node for the current node
    std::transform(tp_gf.begin(), tp_gf.end(), tp_gf_at_node.begin(),
                   [&id](const T& vec) { return vec[id]; });
    const auto loc_temperature = tp_gf_at_node[0];
    const auto& epsilon = 1.e-10;
    int i = 0;
    int j = 0;
    int k = 0;
    // const auto& elem = std::get<0>(chemical_system[0]);
    const auto x1 = tp_gf_at_node[2];
    const auto x2 = tp_gf_at_node[4];

    auto lower_x1 = std::lower_bound(tab_x1.begin(), tab_x1.end(), x1);

    if (lower_x1 != tab_x1.begin() && (lower_x1 + 1) != tab_x1.end()) {
      j = lower_x1 - tab_x1.begin() - 1;
    } else if ((lower_x1 + 1) == tab_x1.end()) {
      j = tab_x1.size() - 1;
    }

    auto lower_x2 = std::lower_bound(tab_x2.begin(), tab_x2.end(), x2);

    if (lower_x2 != tab_x2.begin() && (lower_x2 + 1) != tab_x2.end()) {
      k = lower_x2 - tab_x2.begin() - 1;
    } else if ((lower_x2 + 1) == tab_x2.end()) {
      k = tab_x2.size() - 1;
    }

    auto lower_t = std::lower_bound(tab_T.begin(), tab_T.end(), loc_temperature);

    if (lower_t != tab_T.begin() && (lower_t + 1) != tab_T.end()) {
      i = lower_t - tab_T.begin() - 1;
    } else if ((lower_t + 1) == tab_T.end()) {
      i = tab_T.size() - 1;
    }

    double t = (loc_temperature - tab_T[i]) / (tab_T[i + 1] - tab_T[i]);
    double u = (x1 - tab_x1[j]) / (tab_x1[j + 1] - tab_x1[j]);
    double v = (x2 - tab_x2[k]) / (tab_x2[k + 1] - tab_x2[k]);

    // Energy
    for (const auto& energy_name : energy_names) {
      const auto& key = std::make_tuple(id, phase, energy_name);

      this->energies_of_phases_[key] =
          (1 - t) * (1 - u) * (1 - v) * array_g[i][j][k] +
          t * (1 - u) * (1 - v) * array_g[i + 1][j][k] +
          (1 - t) * u * (1 - v) * array_g[i][j + 1][k] +
          (1 - t) * (1 - u) * v * array_g[i][j][k + 1] +
          t * u * (1 - v) * array_g[i + 1][j + 1][k] + t * (1 - u) * v * array_g[i + 1][j][k + 1] +
          (1 - t) * u * v * array_g[i][j + 1][k + 1] + t * u * v * array_g[i + 1][j + 1][k + 1];
    }

    for (std::size_t id_elem = 0; id_elem < chemical_system.size(); ++id_elem) {
      const auto& elem = std::get<0>(chemical_system[id_elem]);
      std::cout << id_elem << "  " <<  elem << std::endl;
      // Chemical potentials
      this->chemical_potentials_[std::make_tuple(id, elem)] =
          (1 - t) * (1 - u) * (1 - v) * array_mu[elem][i][j][k] +
          t * (1 - u) * (1 - v) * array_mu[elem][i + 1][j][k] +
          (1 - t) * u * (1 - v) * array_mu[elem][i][j + 1][k] +
          (1 - t) * (1 - u) * v * array_mu[elem][i][j][k + 1] +
          t * u * (1 - v) * array_mu[elem][i + 1][j + 1][k] +
          t * (1 - u) * v * array_mu[elem][i + 1][j][k + 1] +
          (1 - t) * u * v * array_mu[elem][i][j + 1][k + 1] +
          t * u * v * array_mu[elem][i + 1][j + 1][k + 1];
      // Molar fraction
      const auto& key = std::make_tuple(id, phase, elem);

      if (elem == "O") {
        this->elem_mole_fraction_by_phase_[key] = x1;
      }
      if (elem == "U") {
        this->elem_mole_fraction_by_phase_[key] = x2;
      }

      // Mobilities
      this->mobilities_[key] =
          (1 - t) * (1 - u) * (1 - v) * std::exp(array_mobility[elem][i][j][k]) +
          t * (1 - u) * (1 - v) * std::exp(array_mobility[elem][i + 1][j][k]) +
          (1 - t) * u * (1 - v) * std::exp(array_mobility[elem][i][j + 1][k]) +
          (1 - t) * (1 - u) * v * std::exp(array_mobility[elem][i][j][k + 1]) +
          t * u * (1 - v) * std::exp(array_mobility[elem][i + 1][j + 1][k]) +
          t * (1 - u) * v * std::exp(array_mobility[elem][i + 1][j][k + 1]) +
          (1 - t) * u * v * std::exp(array_mobility[elem][i][j + 1][k + 1]) +
          t * u * v * std::exp(array_mobility[elem][i + 1][j + 1][k + 1]);
    }
  }
}
/**
 * @brief Check the consistency of outputs required for the current Calphad problem
 *
 * @tparam T
 * @param output_system
 */
template <typename T>
void GeneralMultiParamsTabulation<T>::check_variables_consistency(
    std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system) {
  for (auto& [output_infos, output_value] : output_system) {
    const std::string& output_type = output_infos.back();

    // Fill output with the relevant values
    switch (calphad_outputs::from(output_type)) {
      case calphad_outputs::gm:
      case calphad_outputs::h:
      case calphad_outputs::hm: {
        MFEM_VERIFY(false, "GeneralMultiParamsTabulation is only built for mu, x and g.\n");

        SlothInfo::debug("Output not available for this Calphad problem: ", output_type);
        break;
      }
    }
  }
}
/**
 * @brief Finalization actions (free memory)
 *
 * @tparam T
 */
template <typename T>
void GeneralMultiParamsTabulation<T>::finalize() {}

/**
 * @brief Destroy the Binary Melting< T>:: Binary Melting object
 *
 * @tparam T
 */
template <typename T>
GeneralMultiParamsTabulation<T>::~GeneralMultiParamsTabulation() {}

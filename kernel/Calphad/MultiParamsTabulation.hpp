/**
 * @file MultiParamsTabulation.hpp
 * @author cp273896  (clement.plumecocq@cea.fr)
 * @brief Tabulated thermodynamic description for an binary isothermal system
 * @version 0.1
 * @date 2025-03-17
 *
 * Copyright CEA (c) 2025
 *
 */


#include <algorithm>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <boost/multi_array.hpp>

#include <H5Cpp.h>


#include "Calphad/CalphadBase.hpp"
#include "Calphad/CalphadUtils.hpp"
#include "Inputs/HDF54Sloth.hpp"
#include "Interpolators/MultiLinearInterpolator.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#pragma once

template <typename T>
class MultiParamsTabulation : CalphadBase<T> {
 private:
  std::vector<double> tab_T;
  std::vector<double> tab_x;
  std::map<std::string, std::vector<double>> tab_parameters;
  boost::multi_array<double, 2> array_g;
  std::map<std::string, boost::multi_array<double, 2>> array_mu;
  std::map<std::string, boost::multi_array<double, 2>> array_mobility;

  inline static std::function<double(double)> scalling_func_mob;
  inline static std::function<double(double)> scalling_func_energy;
  inline static std::function<double(double)> scalling_func_potentials;

  std::unique_ptr<CalphadUtils<T>> CU_;
  void compute(
      const size_t nb_nodes, const std::vector<T>& tp_gf,
      const std::vector<std::tuple<std::string, std::string>>& chemical_system,
      std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system);

  void check_variables_consistency(
      std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system);

 public:
  explicit MultiParamsTabulation(const Parameters& params);

  void initialize() override;

  void execute(const int dt, const std::vector<T>& aux_gf,
               const std::vector<std::tuple<std::string, std::string>>& chemical_system,
               std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>&
                   output_system) override;
  void finalize() override;

  ////////////////////////////////

  void get_parameters() override;

  ////////////////////////////////

  ~MultiParamsTabulation();
};

////////////////////////////////
////////////////////////////////

/**
 * @brief Get the parameters associated with the MultiParamsTabulation object
 *
 * @tparam T
 */
template <typename T>
void MultiParamsTabulation<T>::get_parameters() {
  this->description_ = this->params_.template get_param_value_or_default<std::string>(
      "description",
      "Analytical thermodynamic description for an ideal solution using tabulated data ");
  this->scalling_func_mob =
      this->params_.template get_param_value_or_default<std::function<double(double)>>(
          "scalling_func_mobilities", [](double v) { return std::exp(v); });
  this->scalling_func_energy =
      this->params_.template get_param_value_or_default<std::function<double(double)>>(
          "scalling_func_energy", [](double v) { return v; });
  this->scalling_func_potentials =
      this->params_.template get_param_value_or_default<std::function<double(double)>>(
          "scalling_func_potentials", [](double v) { return v; });
}

////////////////////////////////
////////////////////////////////
/**
 * @brief Construct a new MultiParamsTabulation::MultiParamsTabulation object
 *
 * @param params
 */
template <typename T>
MultiParamsTabulation<T>::MultiParamsTabulation(const Parameters& params) : CalphadBase<T>(params) {
  this->CU_ = std::make_unique<CalphadUtils<T>>();
  this->get_parameters();
}

/**
 * @brief Initialization of the thermodynamic calculation
 *
 * @tparam T
 */
template <typename T>
void MultiParamsTabulation<T>::initialize() {
  Catch_Time_Section("MultiParamsTabulation::initialize");

  std::string filename = this->params_.template get_param_value_or_default<std::string>(
      "data filename", "no input file");

  H5std_string DATASET_NAME;
  const H5std_string FILE_NAME(filename);

  std::vector<std::string> list_of_element =
      this->params_.template get_param_value_or_default<std::vector<std::string>>(
          "list_of_elements", {"O", "U"});
  std::vector<std::string> list_of_dataset_name_mu =
      this->params_.template get_param_value_or_default<std::vector<std::string>>(
          "list_of_dataset_name_mu", {"mu_O", "mu_U"});
  std::vector<std::string> list_of_dataset_name_mob =
      this->params_.template get_param_value_or_default<std::vector<std::string>>(
          "list_of_dataset_name_mob", {"Mo", "Mu"});
  std::vector<std::string> list_of_dataset_name_energies =
      this->params_.template get_param_value_or_default<std::vector<std::string>>(
          "list_of_dataset_name_energies", {"G"});
  std::vector<std::string> list_of_dataset_tabulation_parameters =
      this->params_.template get_param_value_or_default<std::vector<std::string>>(
          "list_of_dataset_tabulation_parameters", {"xO", "T"});

  for (size_t i = 0; i < list_of_element.size(); i++) {
    boost::multi_array<double, 2> myArray;

    DATASET_NAME = list_of_dataset_name_mu[i];
    HDF54Sloth<2>::get_data_from_HDF5(FILE_NAME, DATASET_NAME, myArray);
    array_mu.emplace(list_of_element[i], myArray);

    DATASET_NAME = list_of_dataset_name_mob[i];
    HDF54Sloth<2>::get_data_from_HDF5(FILE_NAME, DATASET_NAME, myArray);
    array_mobility.emplace(list_of_element[i], myArray);
  }

  for (size_t i = 0; i < list_of_dataset_name_energies.size(); i++) {
    DATASET_NAME = list_of_dataset_name_energies[i];
    HDF54Sloth<2>::get_data_from_HDF5(FILE_NAME, DATASET_NAME, array_g);
  }

  // for (size_t i = 0; i < list_of_dataset_tabulation_parameters.size(); i++) {
  //   std::vector<double> temp_vec;
  //   DATASET_NAME = list_of_dataset_tabulation_parameters[i];
  //   HDF54Sloth<1>::get_data_from_HDF5(FILE_NAME, DATASET_NAME, temp_vec);
  //   tab_parameters.emplace(list_of_dataset_tabulation_parameters[i], temp_vec);
  // }
  DATASET_NAME = "xO";
  HDF54Sloth<1>::get_data_from_HDF5(FILE_NAME, DATASET_NAME, tab_x);
  DATASET_NAME = "T";
  HDF54Sloth<1>::get_data_from_HDF5(FILE_NAME, DATASET_NAME, tab_T);
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
void MultiParamsTabulation<T>::execute(
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
void MultiParamsTabulation<T>::compute(
    const size_t nb_nodes, const std::vector<T>& tp_gf,
    const std::vector<std::tuple<std::string, std::string>>& chemical_system,
    std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system) {
  Catch_Time_Section("MultiParamsTabulation::compute");

  // Let us assume an ideal mixing solution
  const std::string& phase = "C1_MO2";
  const std::vector<std::string> energy_names = {"G"};
  std::vector<double> tp_gf_at_node(tp_gf.size());

  std::vector<int> sorted_n_t_p = this->CU_->sort_nodes(tp_gf[0], tp_gf[1], "No", "No");
  const std::size_t N = 2;

  // Process CALPHAD calculations for each node
  for (const auto& id : sorted_n_t_p) {
    // Populate tp_gf_at_node for the current node
    std::transform(tp_gf.begin(), tp_gf.end(), tp_gf_at_node.begin(),
                   [&id](const T& vec) { return vec[id]; });
    const auto loc_temperature = tp_gf_at_node[0];
    const auto& epsilon = 1.e-10;

    const auto x = tp_gf_at_node[2];

    std::array<double, N> point_to_interpolate = {tp_gf_at_node[0], tp_gf_at_node[2]};
    std::array<std::vector<double>, N> grid_values = {tab_T, tab_x};
    std::array<std::size_t, N> lower_indices;

    for (size_t i = 0; i < N; i++) {
      std::size_t index;
      auto lower_index =
          std::lower_bound(grid_values[i].begin(), grid_values[i].end(), point_to_interpolate[i]);

      if (lower_index != grid_values[i].begin() && (lower_index + 1) != grid_values[i].end()) {
        index = lower_index - grid_values[i].begin() - 1;
      } else if ((lower_index + 1) == grid_values[i].end()) {
        index = grid_values[i].size() - 1;
      }
      lower_indices[i] = index;
    }

    std::array<double, N> alpha;
    alpha = MultiLinearInterpolator<double, N>::computeInterpolationCoefficients(
        point_to_interpolate, lower_indices, grid_values);

    // Energy
    for (const auto& energy_name : energy_names) {
      const auto& key = std::make_tuple(id, phase, energy_name);

      this->energies_of_phases_[key] = MultiLinearInterpolator<double, N>::computeInterpolation(
          lower_indices, alpha, array_g, this->scalling_func_energy);
    }

    for (std::size_t id_elem = 0; id_elem < chemical_system.size(); ++id_elem) {
      const auto& elem = std::get<0>(chemical_system[id_elem]);
      // Chemical potentials
      this->chemical_potentials_[std::make_tuple(id, elem)] =
          MultiLinearInterpolator<double, N>::computeInterpolation(
              lower_indices, alpha, array_mu[elem], this->scalling_func_potentials);

      // Molar fraction
      const auto& key = std::make_tuple(id, phase, elem);
      this->elem_mole_fraction_by_phase_[key] = x;
      this->mobilities_[key] = MultiLinearInterpolator<double, N>::computeInterpolation(
          lower_indices, alpha, array_mobility[elem], this->scalling_func_mob);
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
void MultiParamsTabulation<T>::check_variables_consistency(
    std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system) {
  for (auto& [output_infos, output_value] : output_system) {
    const std::string& output_type = output_infos.back();

    // Fill output with the relevant values
    switch (calphad_outputs::from(output_type)) {
      case calphad_outputs::gm:
      case calphad_outputs::h:
      case calphad_outputs::hm: {
        MFEM_VERIFY(false, "MultiParamsTabulation is only built for mu, x and g.\n");

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
void MultiParamsTabulation<T>::finalize() {}

/**
 * @brief Destroy the Binary Melting< T>:: Binary Melting object
 *
 * @tparam T
 */
template <typename T>
MultiParamsTabulation<T>::~MultiParamsTabulation() {}

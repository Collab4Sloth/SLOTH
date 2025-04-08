/**
 * @file MultiParamsTabulation.hpp
 * @author cp273896 (clement.plumecocq@cea.fr)
 * @brief Tabulated thermodynamic description for an binary isothermal system
 * @version 0.1
 * @date 2025-03-17
 *
 * Copyright CEA (c)
 * 2025
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
#include "Interpolators/MultiLinearInterpolator.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#pragma once

template <typename T, std::size_t INTERP_DIM = 2>
class MultiParamsTabulation : CalphadBase<T> {
 private:
  std::string filename;

  boost::multi_array<double, INTERP_DIM> array_g;

  std::map<std::string, boost::multi_array<double, INTERP_DIM>> array_mu;
  std::map<std::string, boost::multi_array<double, INTERP_DIM>> array_mobility;

  std::array<std::vector<double>, INTERP_DIM> grid_values;

  std::vector<std::size_t> list_of_aux_gf_index_for_tabulation;

  std::vector<std::string> list_of_dataset_tabulation_parameters;
  std::vector<std::string> list_of_dataset_name_energies;
  std::vector<std::string> list_of_dataset_name_mob;
  std::vector<std::string> list_of_dataset_name_mu;
  std::vector<std::string> list_of_element;

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
 * @tparam T
 */
template <typename T, std::size_t INTERP_DIM>
void MultiParamsTabulation<T, INTERP_DIM>::get_parameters() {
  this->description_ = this->params_.template get_param_value_or_default<std::string>(
      "description",
      "Analytical thermodynamic description for an ideal solution using tabulated data");
  this->scalling_func_mob =
      this->params_.template get_param_value_or_default<std::function<double(double)>>(
          "scalling_func_mobilities", [](double v) { return std::exp(v); });
  this->scalling_func_energy =
      this->params_.template get_param_value_or_default<std::function<double(double)>>(
          "scalling_func_energy", [](double v) { return v; });
  this->scalling_func_potentials =
      this->params_.template get_param_value_or_default<std::function<double(double)>>(
          "scalling_func_potentials", [](double v) { return v; });
  this->list_of_element =
      this->params_.template get_param_value_or_default<std::vector<std::string>>(
          "list_of_elements", {"O", "U"});
  this->list_of_dataset_name_mu =
      this->params_.template get_param_value_or_default<std::vector<std::string>>(
          "list_of_dataset_name_mu", {"mu_O", "mu_U"});
  this->list_of_dataset_name_mob =
      this->params_.template get_param_value_or_default<std::vector<std::string>>(
          "list_of_dataset_name_mob", {"Mo", "Mu"});
  this->list_of_dataset_name_energies =
      this->params_.template get_param_value_or_default<std::vector<std::string>>(
          "list_of_dataset_name_energies", {"G"});
  this->list_of_dataset_tabulation_parameters =
      this->params_.template get_param_value_or_default<std::vector<std::string>>(
          "list_of_dataset_tabulation_parameters", {"T", "xO"});
  this->list_of_aux_gf_index_for_tabulation =
      this->params_.template get_param_value_or_default<std::vector<std::size_t>>(
          "list_of_aux_gf_index_for_tabulation", {0, 2});
  this->filename = this->params_.template get_param_value_or_default<std::string>("data filename",
                                                                                  "no input file");
}

////////////////////////////////
////////////////////////////////
/**
 * @brief Construct a new MultiParamsTabulation::MultiParamsTabulation object
 *
 * @param params
 */
template <typename T, std::size_t INTERP_DIM>
MultiParamsTabulation<T, INTERP_DIM>::MultiParamsTabulation(const Parameters& params)
    : CalphadBase<T>(params) {
  this->CU_ = std::make_unique<CalphadUtils<T>>();
  this->get_parameters();
}

/**
 * @brief Initialization of the thermodynamic calculation
 * @tparam T
 */
template <typename T, std::size_t INTERP_DIM>
void MultiParamsTabulation<T, INTERP_DIM>::initialize() {
  Catch_Time_Section(
      "MultiParamsT"
      "abulation::"
      "initialize");

  H5std_string DATASET_NAME;
  const H5std_string FILE_NAME(this->filename);

  for (size_t i = 0; i < list_of_element.size(); i++) {
    boost::multi_array<double, INTERP_DIM> myArray;

    DATASET_NAME = list_of_dataset_name_mu[i];
    HDF54Sloth<INTERP_DIM>::get_data_from_HDF5(FILE_NAME, DATASET_NAME, myArray);
    array_mu.emplace(list_of_element[i], myArray);

    DATASET_NAME = list_of_dataset_name_mob[i];
    HDF54Sloth<INTERP_DIM>::get_data_from_HDF5(FILE_NAME, DATASET_NAME, myArray);
    array_mobility.emplace(list_of_element[i], myArray);
  }

  for (size_t i = 0; i < list_of_dataset_name_energies.size(); i++) {
    DATASET_NAME = list_of_dataset_name_energies[i];
    HDF54Sloth<INTERP_DIM>::get_data_from_HDF5(FILE_NAME, DATASET_NAME, array_g);
  }

  for (size_t i = 0; i < list_of_dataset_tabulation_parameters.size(); i++) {
    std::vector<double> temp_vec;
    DATASET_NAME = list_of_dataset_tabulation_parameters[i];
    HDF54Sloth<1>::get_data_from_HDF5(FILE_NAME, DATASET_NAME, temp_vec);
    grid_values[i] = temp_vec;
  }
}

/**
 * @brief Main method to calculate equilibrium states
 * @tparam T
 * @param dt
 * @param aux_gf
 * @param
 * chemical_system
 * @param
 * output_system
 */
template <typename T, std::size_t INTERP_DIM>
void MultiParamsTabulation<T, INTERP_DIM>::execute(
    const int dt, const std::vector<T>& tp_gf,
    const std::vector<std::tuple<std::string, std::string>>& chemical_system,
    std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system) {
  const size_t nb_nodes = this->CU_->get_size(tp_gf[0]);
  this->clear_containers();
  if (dt == 1) {
    this->check_variables_consistency(output_system);
  }

  this->compute(nb_nodes, tp_gf, chemical_system, output_system);

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
template <typename T, std::size_t INTERP_DIM>
void MultiParamsTabulation<T, INTERP_DIM>::compute(
    const size_t nb_nodes, const std::vector<T>& tp_gf,
    const std::vector<std::tuple<std::string, std::string>>& chemical_system,
    std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system) {
  Catch_Time_Section(
      "MultiParamsT"
      "abulation::"
      "compute");

  const std::string& phase = "C1_MO2";
  const std::vector<std::string> energy_names = {"G"};
  std::vector<double> tp_gf_at_node(tp_gf.size());

  std::vector<int> sorted_n_t_p = this->CU_->sort_nodes(tp_gf[0], tp_gf[1], "No", "No");
  const std::size_t N = INTERP_DIM;

  std::array<double, N> point_to_interpolate;
  std::array<std::size_t, N> lower_indices;
  std::array<double, N> alpha;

  for (const auto& id : sorted_n_t_p) {
    std::transform(tp_gf.begin(), tp_gf.end(), tp_gf_at_node.begin(),
                   [&id](const T& vec) { return vec[id]; });

    for (size_t i = 0; i < N; i++) {
      point_to_interpolate[i] = tp_gf_at_node[list_of_aux_gf_index_for_tabulation[i]];
    }

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

    alpha = MultiLinearInterpolator<double, N>::computeInterpolationCoefficients(
        point_to_interpolate, lower_indices, grid_values);

    // Energies
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

      // Molar
      // fractions
      const auto& key = std::make_tuple(id, phase, elem);
      this->elem_mole_fraction_by_phase_[key] = tp_gf_at_node[id_elem + 2];

      // Mobilities
      this->mobilities_[key] = MultiLinearInterpolator<double, N>::computeInterpolation(
          lower_indices, alpha, array_mobility[elem], this->scalling_func_mob);
    }
  }
}
/**
 * @brief Check the consistency of outputs required for the current Calphad problem
 * @tparam T
 * @param output_system
 */
template <typename T, std::size_t INTERP_DIM>
void MultiParamsTabulation<T, INTERP_DIM>::check_variables_consistency(
    std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system) {
  for (auto& [output_infos, output_value] : output_system) {
    const std::string& output_type = output_infos.back();

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
template <typename T, std::size_t INTERP_DIM>
void MultiParamsTabulation<T, INTERP_DIM>::finalize() {}

/**
 * @brief Destroy the Binary Melting< T>::Binary Melting object
 * @tparam T
 */
template <typename T, std::size_t INTERP_DIM>
MultiParamsTabulation<T, INTERP_DIM>::~MultiParamsTabulation() {}

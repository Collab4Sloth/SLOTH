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
#include "kernel/Utils/Utils.hpp"
#pragma once

template <typename T>
class MultiParamsTabulation : public CalphadBase<T> {
 private:
  std::map<std::string, unsigned int> idx_elem_;
  std::vector<std::string> input_composition_order_;
  std::size_t INTERP_DIM;
  std::string filename;

  FlattenedTensor<double> flatten_g;
  std::map<std::string, FlattenedTensor<double>> flatten_mu;
  std::map<std::string, FlattenedTensor<double>> flatten_mobility;

  std::vector<std::vector<double>> grid_values;

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
  void compute(const size_t nb_nodes, const std::vector<T>& tp_gf,
               const std::vector<std::tuple<std::string, std::string>>& chemical_system);

  void check_variables_consistency(
      std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system);

 public:
  constexpr explicit MultiParamsTabulation(const Parameters& params);

  constexpr MultiParamsTabulation(const Parameters& params, bool is_KKS);

  void initialize(
      const std::vector<std::tuple<std::string, std::string>>& sorted_chemical_system) override;

  void execute(const int dt, const std::set<int>& list_nodes, const std::vector<T>& tp_gf,
               const std::vector<std::tuple<std::string, std::string>>& chemical_system,
               std::optional<std::vector<std::tuple<std::string, std::string, double>>>
                   status_phase = std::nullopt) override;
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
template <typename T>
void MultiParamsTabulation<T>::get_parameters() {
  // this->description_ = this->params_.template get_param_value_or_default<std::string>(
  //     "description",
  //     "Analytical thermodynamic description for an ideal solution using tabulated data");
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
  this->INTERP_DIM = static_cast<std::size_t>(
      this->params_.template get_param_value<int>("dimension_of_interpolation"));

  bool check_parameters_mob_mu =
      (list_of_dataset_name_mob.size() == list_of_dataset_name_mu.size()) &&
      (list_of_dataset_name_mob.size() == list_of_element.size());
  bool check_parameters_tab_params =
      (list_of_aux_gf_index_for_tabulation.size() ==
       list_of_dataset_tabulation_parameters.size()) &&
      (list_of_aux_gf_index_for_tabulation.size() == this->INTERP_DIM);
  MFEM_VERIFY(
      check_parameters_tab_params,
      " Parameters : (1) list_of_aux_gf_index_for_tabulation , (2) "
      "list_of_dataset_tabulation_parameters , (3) dimension_of_interpolation are not consistent");
  MFEM_VERIFY(check_parameters_mob_mu,
              " Parameters : (1) list_of_dataset_name_mob , (2) "
              "list_of_dataset_name_mu , (3) list_of_element are not consistent");
}

////////////////////////////////
////////////////////////////////
/**
 * @brief Construct a new MultiParamsTabulation::MultiParamsTabulation object
 *
 * @param params
 */
template <typename T>
constexpr MultiParamsTabulation<T>::MultiParamsTabulation(const Parameters& params)
    : CalphadBase<T>(params, false) {
  this->CU_ = std::make_unique<CalphadUtils<T>>();
  this->get_parameters();
}

/**
 * @brief Construct a new MultiParamsTabulation::MultiParamsTabulation object
 *
 * @param params
 */
template <typename T>
constexpr MultiParamsTabulation<T>::MultiParamsTabulation(const Parameters& params, bool is_KKS)
    : CalphadBase<T>(params, is_KKS) {
  this->CU_ = std::make_unique<CalphadUtils<T>>();
  this->get_parameters();
}

/**
 * @brief Initialization of the thermodynamic calculation
 * @tparam T
 */
template <typename T>
void MultiParamsTabulation<T>::initialize(
    const std::vector<std::tuple<std::string, std::string>>& sorted_chemical_system) {
  Catch_Time_Section("MultiParamsTabulation::initialize");

  H5std_string DATASET_NAME;
  const H5std_string FILE_NAME(this->filename);

  for (size_t i = 0; i < list_of_element.size(); i++) {
    FlattenedTensor<double> myVec;

    DATASET_NAME = list_of_dataset_name_mu[i];
    HDF54Sloth<std::monostate>::get_data_from_HDF5(FILE_NAME, DATASET_NAME, myVec);
    myVec.apply_scalling(this->scalling_func_potentials);
    flatten_mu.emplace(list_of_element[i], myVec);

    DATASET_NAME = list_of_dataset_name_mob[i];
    HDF54Sloth<std::monostate>::get_data_from_HDF5(FILE_NAME, DATASET_NAME, myVec);
    myVec.apply_scalling(this->scalling_func_mob);
    flatten_mobility.emplace(list_of_element[i], myVec);
  }

  for (size_t i = 0; i < list_of_dataset_name_energies.size(); i++) {
    DATASET_NAME = list_of_dataset_name_energies[i];
    HDF54Sloth<std::monostate>::get_data_from_HDF5(FILE_NAME, DATASET_NAME, flatten_g);
    flatten_g.apply_scalling(this->scalling_func_energy);
  }
  grid_values.resize(this->INTERP_DIM);
  for (size_t i = 0; i < list_of_dataset_tabulation_parameters.size(); i++) {
    std::vector<double> temp_vec;
    DATASET_NAME = list_of_dataset_tabulation_parameters[i];
    HDF54Sloth<std::monostate>::get_data_from_HDF5(FILE_NAME, DATASET_NAME, temp_vec);
    grid_values[i] = temp_vec;
  }

  // Check input composition consistency (cf CalphadInformedNeuralNetwork.hpp)
  if (this->params_.template has_parameter("InputCompositionOrder")) {
    this->input_composition_order_ =
        this->params_.template get_param_value<vString>("InputCompositionOrder");
  }

  for (int i = 0; i < sorted_chemical_system.size(); i++) {
    const std::string& elem = std::get<0>(sorted_chemical_system[i]);
    this->idx_elem_.emplace(elem, i);
  }
  if (!this->input_composition_order_.empty()) {
    std::vector<std::string> v_elem;
    for (const auto& tup : sorted_chemical_system) {
      const std::string& elem = std::get<0>(tup);
      v_elem.emplace_back(elem);
      this->idx_elem_[elem] = std::distance(this->input_composition_order_.begin(),
                                            std::find(this->input_composition_order_.begin(),
                                                      this->input_composition_order_.end(), elem));
    }
    // Sort before comparison
    std::ranges::sort(v_elem);
    auto sort_input_composition_order = this->input_composition_order_;
    std::ranges::sort(sort_input_composition_order);

    MFEM_VERIFY(std::ranges::equal(v_elem, sort_input_composition_order),
                "Error: InputCompositionOrder is not consistent with composition deduced from "
                "auxiliary variables. Please check your data");
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
template <typename T>
void MultiParamsTabulation<T>::execute(
    const int dt, const std::set<int>& list_nodes, const std::vector<T>& tp_gf,
    const std::vector<std::tuple<std::string, std::string>>& chemical_system,
    std::optional<std::vector<std::tuple<std::string, std::string, double>>> status_phase) {
  // const size_t nb_nodes = this->CU_->get_size(tp_gf[0]);
  const size_t nb_nodes = list_nodes.size();
  this->clear_containers();
  // if (dt == 1) {
  //   this->check_variables_consistency(output_system);
  // }

  this->compute(nb_nodes, tp_gf, chemical_system);

  // this->update_outputs(nb_nodes, output_system);
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
    const std::vector<std::tuple<std::string, std::string>>& chemical_system) {
  Catch_Time_Section("MultiParamsTabulation::compute");

  const std::string& phase = "C1_MO2";
  const std::vector<std::string> energy_names = {"G"};
  std::vector<double> tp_gf_at_node(tp_gf.size());

  std::vector<int> sorted_n_t_p = this->CU_->sort_nodes(tp_gf[0], tp_gf[1], "No", "No");
  std::size_t N = this->INTERP_DIM;

  std::vector<double> point_to_interpolate(N);
  std::vector<std::size_t> lower_indices(N);
  std::vector<double> alpha(N);
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

    alpha = MultiLinearInterpolator<double>::computeInterpolationCoefficients(
        N, point_to_interpolate, lower_indices, grid_values);

    // Energies
    for (const auto& energy_name : energy_names) {
      const auto& key = std::make_tuple(id, phase, energy_name);

      this->energies_of_phases_[key] =
          MultiLinearInterpolator<double>::computeInterpolation(N, lower_indices, alpha, flatten_g);
    }

    for (std::size_t id_elem = 0; id_elem < chemical_system.size(); ++id_elem) {
      const auto& elem = std::get<0>(chemical_system[id_elem]);
      // Chemical potentials

      this->chemical_potentials_[std::make_tuple(id, elem)] =
          MultiLinearInterpolator<double>::computeInterpolation(N, lower_indices, alpha,
                                                                flatten_mu[elem]);

      // Molar fractions
      const auto& key = std::make_tuple(id, phase, elem);
      this->elem_mole_fraction_by_phase_[key] = tp_gf_at_node[this->idx_elem_[elem] + 2];

      // Mobilities
      this->mobilities_[key] = MultiLinearInterpolator<double>::computeInterpolation(
          N, lower_indices, alpha, flatten_mobility[elem]);
    }
  }
}
/**
 * @brief Check the consistency of outputs required for the current Calphad problem
 * @tparam T
 * @param output_system
 */
template <typename T>
void MultiParamsTabulation<T>::check_variables_consistency(
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
template <typename T>
void MultiParamsTabulation<T>::finalize() {}

/**
 * @brief Destroy the Binary Melting< T>::Binary Melting object
 * @tparam T
 */
template <typename T>
MultiParamsTabulation<T>::~MultiParamsTabulation() {}

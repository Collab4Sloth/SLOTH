/**
 * @file CalphadInformedNeuralNetwork.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Calphad-informed neuronal network used to predict Calphad contributions
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
#include <torch/script.h>
#include <torch/torch.h>

#include <algorithm>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "Calphad/CalphadBase.hpp"
#include "Calphad/CalphadUtils.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"

#pragma once

/**
 * @brief Calphad calculation performed by neural network
 * Dedicated to simulations with multicomponent inter-diffusion
 * @tparam T
 */
template <typename T>
class CalphadInformedNeuralNetwork : public CalphadBase<T> {
 private:
  // Name of a phase for calphad calculations in only one phase
  std::string given_phase_;
  // List of neural network by phase and associated containers : chemical potentials
  std::unordered_map<std::string, std::unique_ptr<torch::jit::script::Module>>
      chemical_potentials_nn_;
  std::unordered_map<std::string, std::string> chemical_potentials_neural_network_;
  std::unordered_map<std::string, int> index_chemical_potentials_neural_network_;
  // List of neural network by phase and associated containers  : mobilities
  std::unordered_map<std::string, std::unique_ptr<torch::jit::script::Module>> mobilities_nn_;
  std::unordered_map<std::string, std::string> mobilities_neural_network_;
  std::unordered_map<std::string, int> index_mobilities_neural_network_;

  // List of neural network by phase and associated containers  : energies (G, GM, H, HM)
  std::unordered_map<std::string, std::unique_ptr<torch::jit::script::Module>> energies_nn_;
  std::unordered_map<std::string, std::string> energies_neural_network_;
  std::unordered_map<std::string, int> index_energies_neural_network_;

  // List of flag to identify the calphad contributions calculated by neural network
  bool has_nn_potentials_{false};
  bool has_nn_mobilities_{false};
  bool has_nn_energies_{false};
  // List of flag to know if mobilities and/or energies are calculated by the neural network of the
  // chemical potentials or their own meta-model
  bool own_mobilities_model_{false};
  bool own_energies_model_{false};

  // Flag to identify if the pressure is an input of neural networks (rule is the same for all
  // meta-models)
  bool model_built_with_pressure_{false};

  // Objects used to manage the composition and energies in all meta-models
  // Name of the element removed from inputs in all meta-models
  bool element_removed_from_inputs_{false};
  std::string element_removed_from_nn_inputs_;
  // Order of the chemical elements in the inputs and outputs for all meta-models
  std::vector<std::string> input_composition_order_;
  std::map<std::string, unsigned int> idx_elem_;
  // Factor used in case of inputs in moles (instead of molar fractions)
  double input_composition_factor_;
  // Order of the energies in output of the meta-model
  std::vector<std::string> input_energies_order_;

  std::unique_ptr<CalphadUtils<T>> CU_;
  void compute(const std::set<int>& list_nodes, const std::vector<T>& tp_gf,
               const std::vector<std::tuple<std::string, std::string>>& chemical_system,
               const std::string& given_phase);

  void check_variables_consistency(
      std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system);

 public:
  constexpr explicit CalphadInformedNeuralNetwork(const Parameters& params);
  constexpr CalphadInformedNeuralNetwork(const Parameters& params, bool is_KKS);

  void initialize(
      const std::vector<std::tuple<std::string, std::string>>& sorted_chemical_system) override;

  void execute(const int dt, const std::set<int>& list_nodes, const std::vector<T>& aux_gf,
               const std::vector<std::tuple<std::string, std::string>>& chemical_system,
               std::optional<std::vector<std::tuple<std::string, std::string, double>>>
                   status_phase = std::nullopt) override;

  void finalize() override;

  ////////////////////////////////

  void get_parameters() override;

  ////////////////////////////////

  ~CalphadInformedNeuralNetwork();
};

////////////////////////////////
////////////////////////////////

/**
 * @brief Get the parameters associated with the CalphadInformedNeuralNetwork object
 *
 * @tparam T
 */
template <typename T>
void CalphadInformedNeuralNetwork<T>::get_parameters() {
  CalphadBase<T>::get_parameters();

  // Model and Index for Chemical potentials
  if (this->params_.has_parameter("ChemicalPotentialsNeuralNetwork")) {
    this->has_nn_potentials_ = true;
    std::set<std::string> set_mu_phase;
    const auto& vChemicalPotentialNeuralNetwork =
        this->params_.template get_param_value<vTupleStringString>(
            "ChemicalPotentialsNeuralNetwork");
    for (const auto& v : vChemicalPotentialNeuralNetwork) {
      const auto& [mu_model, mu_phase] = v;
      check_file_existence(mu_model);
      this->chemical_potentials_neural_network_.emplace(mu_phase, mu_model);
      set_mu_phase.insert(mu_phase);
    }

    if (this->params_.has_parameter("ChemicalPotentialsNeuralNetworkIndex")) {
      const auto& vChemicalPotentialNeuralNetworkIndex =
          this->params_.template get_param_value<vTupleStringInt>(
              "ChemicalPotentialsNeuralNetworkIndex");
      for (const auto& v : vChemicalPotentialNeuralNetworkIndex) {
        const auto& [mu_phase, index_mu_model] = v;
        if (set_mu_phase.find(mu_phase) != set_mu_phase.end()) {
          this->index_chemical_potentials_neural_network_.emplace(mu_phase, index_mu_model);
        } else {
          const std::string& error_msg =
              "Error while loading index of the neuronal network for chemical potentials. Please "
              "check the name of the phases in the parameters ChemicalPotentialsNeuralNetwork and "
              "ChemicalPotentialsNeuralNetworkIndex.";
          mfem::mfem_error(error_msg.c_str());
        }
      }

    } else {
      for (const auto& phase : set_mu_phase) {
        this->index_chemical_potentials_neural_network_.emplace(phase, 0);
      }
    }
    set_mu_phase.clear();
  }

  // Model and Index for Mobilities
  if (this->params_.has_parameter("MobilitiesNeuralNetwork")) {
    this->has_nn_mobilities_ = true;
    std::set<std::string> set_mob_phase;
    const auto& vMobilitiesNeuralNetwork =
        this->params_.template get_param_value<vTupleStringString>("MobilitiesNeuralNetwork");
    for (const auto& v : vMobilitiesNeuralNetwork) {
      const auto& [mob_model, mob_phase] = v;
      check_file_existence(mob_model);
      this->mobilities_neural_network_.emplace(mob_phase, mob_model);
      set_mob_phase.insert(mob_phase);
    }

    if (this->params_.has_parameter("MobilitiesNeuralNetworkIndex")) {
      const auto& vMobilitiesNeuralNetworkIndex =
          this->params_.template get_param_value<vTupleStringInt>("MobilitiesNeuralNetworkIndex");
      for (const auto& v : vMobilitiesNeuralNetworkIndex) {
        const auto& [mob_phase, index_mob_model] = v;
        if (set_mob_phase.find(mob_phase) != set_mob_phase.end()) {
          this->index_mobilities_neural_network_.emplace(mob_phase, index_mob_model);
        } else {
          const std::string& error_msg =
              "Error while loading index of the neuronal network for mobilities. Please check the "
              "name of the phases in the parameters MobilitiesNeuralNetwork and "
              "MobilitiesNeuralNetworkIndex.";
          mfem::mfem_error(error_msg.c_str());
        }
      }

    } else {
      for (const auto& phase : set_mob_phase) {
        this->index_mobilities_neural_network_.emplace(phase, 0);
      }
    }
    set_mob_phase.clear();
  }

  // Model and Index Energies
  if (this->params_.has_parameter("EnergiesNeuralNetwork")) {
    this->has_nn_energies_ = true;
    std::set<std::string> set_energies_phase;
    const auto& vEnergiesNeuralNetwork =
        this->params_.template get_param_value<vTupleStringString>("EnergiesNeuralNetwork");
    for (const auto& v : vEnergiesNeuralNetwork) {
      const auto& [nrj_model, nrj_phase] = v;
      check_file_existence(nrj_model);
      this->energies_neural_network_.emplace(nrj_phase, nrj_model);
      set_energies_phase.insert(nrj_phase);
    }

    if (this->params_.has_parameter("EnergiesNeuralNetworkIndex")) {
      const auto& vEnergiesNeuralNetworkIndex =
          this->params_.template get_param_value<vTupleStringInt>("EnergiesNeuralNetworkIndex");
      for (const auto& v : vEnergiesNeuralNetworkIndex) {
        const auto& [nrj_phase, index_nrj_model] = v;
        if (set_energies_phase.find(nrj_phase) != set_energies_phase.end()) {
          this->index_energies_neural_network_.emplace(nrj_phase, index_nrj_model);
        } else {
          const std::string& error_msg =
              "Error while loading index of the neuronal network for energies. Please check the "
              "name of the phases in the parameters EnergiesNeuralNetwork and "
              "MEnergiessNeuralNetworkIndex.";
          mfem::mfem_error(error_msg.c_str());
        }
      }
    } else {
      for (const auto& phase : set_energies_phase) {
        this->index_energies_neural_network_.emplace(phase, 0);
      }
    }
    set_energies_phase.clear();
    this->input_energies_order_ =
        this->params_.template get_param_value<vString>("InputEnergiesOrder");
  }

  // Composition order if different from alphabetical order
  if (this->params_.has_parameter("InputCompositionOrder")) {
    this->input_composition_order_ =
        this->params_.template get_param_value<vString>("InputCompositionOrder");
  }
  this->input_composition_factor_ =
      this->params_.template get_param_value_or_default<double>("InputCompositionFactor", 1.);

  // Check if pressure is an input of the model
  if (this->params_.has_parameter("ModelBuiltWithPressure")) {
    this->model_built_with_pressure_ =
        this->params_.template get_param_value<bool>("ModelBuiltWithPressure");
  }

  // Check if a specific model is used for mobilities and energies.
  // Otherwise the model for chemical potential is used
  if (this->params_.has_parameter("OwnMobilityModel")) {
    this->own_mobilities_model_ = this->params_.template get_param_value<bool>("OwnMobilityModel");
  }
  if (this->params_.has_parameter("OwnEnergyModel")) {
    this->own_energies_model_ = this->params_.template get_param_value<bool>("OwnEnergyModel");
  }

  // Given phase : parameter used for single phase calculation (without KKS module)
  this->given_phase_ =
      this->params_.template get_param_value_or_default<std::string>("GivenPhase", "");

  if (this->params_.has_parameter("element_removed_from_nn_inputs")) {
    this->element_removed_from_nn_inputs_ =
        this->params_.template get_param_value<std::string>("element_removed_from_nn_inputs");
    this->element_removed_from_inputs_ = true;
  }
}

////////////////////////////////
////////////////////////////////
/**
 * @brief Construct a new CalphadInformedNeuralNetwork::CalphadInformedNeuralNetwork object
 *
 * @param params The parameter of the Calphad problem
 */
template <typename T>
constexpr CalphadInformedNeuralNetwork<T>::CalphadInformedNeuralNetwork(const Parameters& params)
    : CalphadBase<T>(params, false) {
  this->CU_ = std::make_unique<CalphadUtils<T>>();

  // torch::set_num_threads(8);

  this->get_parameters();
}

/**
 * @brief Construct a new CalphadInformedNeuralNetwork::CalphadInformedNeuralNetwork object
 *
 * @param params The parameter of the Calphad problem
 * @param is_KKS A boolean flag indicating whether to initialize KKS.
 */
template <typename T>
constexpr CalphadInformedNeuralNetwork<T>::CalphadInformedNeuralNetwork(const Parameters& params,
                                                                        bool is_KKS)
    : CalphadBase<T>(params, is_KKS) {
  this->CU_ = std::make_unique<CalphadUtils<T>>();

  // torch::set_num_threads(8);

  this->get_parameters();
}

/**
 * @brief Initialization of the thermodynamic calculation
 *
 * @tparam T
 */
template <typename T>
void CalphadInformedNeuralNetwork<T>::initialize(
    const std::vector<std::tuple<std::string, std::string>>& sorted_chemical_system) {
  Catch_Time_Section("CalphadInformedNeuralNetwork<T>::initialize");
  at::set_num_interop_threads(1);
  at::set_num_threads(1);
  // torch::AutoGradMode enable_grad(false);
  ////////////////////////////////////////////////////////////
  // Check input composition consistency
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
  ///////////////////////////////////////////////////////////////
  // Deserialize the ScriptModule for chemical potentials
  if (!this->chemical_potentials_neural_network_.empty()) {
    try {
      for (const auto& [phase, mmodel] : this->chemical_potentials_neural_network_) {
        // Deserialize the ScriptModule from a file using torch::jit::load()
        this->chemical_potentials_nn_[phase] =
            std::make_unique<torch::jit::script::Module>(torch::jit::load(mmodel));
      }
    } catch (const c10::Error& e) {
      const std::string& error_msg =
          "Error while loading  neuronal network for chemical potentials. Please  check the "
          "input datafile";
      mfem::mfem_error(error_msg.c_str());
    }
  }

  ///////////////////////////////////////////////////////////////
  // Deserialize the ScriptModules for mobilities by phase
  if (!this->mobilities_neural_network_.empty()) {
    try {
      for (const auto& [phase, mmodel] : this->mobilities_neural_network_) {
        MFEM_VERIFY(phase != "LIQUID", "Error with mobility model deserialize ");

        // Deserialize the ScriptModule from a file using torch::jit::load()
        this->mobilities_nn_[phase] =
            std::make_unique<torch::jit::script::Module>(torch::jit::load(mmodel));
      }
    } catch (const c10::Error& e) {
      const std::string& error_msg =
          "Error while loading  neuronal network for mobilities. Please  check the "
          "input datafile";
      mfem::mfem_error(error_msg.c_str());
    }
  }

  ///////////////////////////////////////////////////////////////
  // Deserialize the ScriptModules for mobilities by phase
  if (!this->energies_neural_network_.empty()) {
    try {
      for (const auto& [phase, mmodel] : this->energies_neural_network_) {
        // Deserialize the ScriptModule from a file using torch::jit::load()
        this->energies_nn_[phase] =
            std::make_unique<torch::jit::script::Module>(torch::jit::load(mmodel));
      }
    } catch (const c10::Error& e) {
      const std::string& error_msg =
          "Error while loading  neuronal network for energies. Please  check the "
          "input datafile";
      mfem::mfem_error(error_msg.c_str());
    }
  }
}

/**
 * @brief Perform equilibrium calculations in a specific phase
 *
 * @tparam T
 * @param dt The current iteration  of the simulation
 * @param list_nodes The list of nodes where calculations are performed
 * @param tp_gf The thermodynamic condition for equilibrium calculations
 * @param chemical_system  The targeted chemical system
 * @param status_phase A vector of phases with their status
 */
template <typename T>
void CalphadInformedNeuralNetwork<T>::execute(
    const int dt, const std::set<int>& list_nodes, const std::vector<T>& tp_gf,
    const std::vector<std::tuple<std::string, std::string>>& chemical_system,
    std::optional<std::vector<std::tuple<std::string, std::string, double>>> status_phase) {
  Catch_Time_Section("CalphadInformedNeuralNetwork<T>::execute");

  // Thermodynamic calculations are done by phase
  std::string phase;
  if (status_phase.has_value()) {
    const auto& phase_vector = *status_phase;
    MFEM_VERIFY(phase_vector.size() == 1,
                "CalphadInformedNeuralNetwork: status phase must contain exactly one element. \n");
    phase = std::get<0>(phase_vector[0]);
  } else {
    if (this->given_phase_.size() > 0) {
      phase = this->given_phase_;
    } else {
      throw std::runtime_error(
          "CalphadInformedNeuralNetwork::execute: GivenPhase parameter must be defined.");
    }
  }

  this->compute(list_nodes, tp_gf, chemical_system, phase);
}

/**
 * @brief Compute the CALPHAD contributions from Neural Networks
 *
 * @tparam T
 * @param tp_gf
 * @param chemical_system
 * @param output_system
 */
template <typename T>
void CalphadInformedNeuralNetwork<T>::compute(
    const std::set<int>& list_nodes, const std::vector<T>& tp_gf,
    const std::vector<std::tuple<std::string, std::string>>& chemical_system,
    const std::string& given_phase) {
  Catch_Time_Section("CalphadInformedNeuralNetwork<T>::compute");

  /////////////////////////////////////////////////////////////////////////////
  // Build tensor of inputs (all nodes)
  /////////////////////////////////////////////////////////////////////////////
  const size_t nb_nodes = list_nodes.size();
  std::vector<double> tp_gf_at_node(tp_gf.size());
  // By default, the number of inputs is the number of chemical element + temperature
  int nb_input = chemical_system.size() + 1;
  // The number of inputs is reduced if the reference element is removed from inputs
  // (consistent when working with molar fractions)
  if (this->element_removed_from_inputs_) {
    nb_input--;
  }
  // The number of inputs is enhanced if the pressure is an input of all meta-models
  if (this->model_built_with_pressure_) {
    nb_input++;
  }

  auto options = torch::TensorOptions().dtype(torch::kDouble);
  auto tensor =
      torch::zeros({static_cast<int64_t>(nb_nodes), static_cast<int64_t>(nb_input)}, options);

  // ia is the index of tensor. ia is not equal to i because the list of nodes is not necessary from
  // 0 to nb_nodes
  int ia = -1;
  for (const auto& i : list_nodes) {
    ia++;
    // Populate tp_gf_at_node for the current node
    std::transform(tp_gf.begin(), tp_gf.end(), tp_gf_at_node.begin(),
                   [&i](const T& vec) { return vec[i]; });

    const double temperature = tp_gf_at_node[0];
    const double pressure = tp_gf_at_node[1];

    tensor[ia][0] = temperature;

    // id_elem starts to 1 if pressure is not an input of the model. Otherwise it starts to 2.
    int id_elem = 1;
    if (this->model_built_with_pressure_) {
      tensor[ia][1] = pressure;
      id_elem = 2;
    }
    // sum is used to calculate the molar fraction of the element removed from the inputs
    double sum_x = 0.0;
    for (std::size_t iel = 0; iel < chemical_system.size(); ++iel) {
      const std::string& elem = std::get<0>(chemical_system[iel]);
      if (this->element_removed_from_inputs_ && elem == this->element_removed_from_nn_inputs_)
        continue;

      int idx_el = this->idx_elem_[elem];
      if ((this->element_removed_from_inputs_) &&
          (this->idx_elem_[elem] > this->idx_elem_[this->element_removed_from_nn_inputs_])) {
        idx_el = this->idx_elem_[elem] - 1;
      }

      tensor[ia][idx_el + id_elem] = tp_gf_at_node[iel + 2] * input_composition_factor_;
      const auto& key = std::make_tuple(i, given_phase, elem);
      this->elem_mole_fraction_by_phase_[key] = tp_gf_at_node[iel + 2];
      if (this->element_removed_from_inputs_) {
        sum_x += this->elem_mole_fraction_by_phase_[key];
      }
    }
    // Molar fraction of the removed element by phase
    if (this->element_removed_from_inputs_) {
      this->elem_mole_fraction_by_phase_[std::make_tuple(
          i, given_phase, this->element_removed_from_nn_inputs_)] = 1.0 - sum_x;
    }
  }

  std::vector<torch::jit::IValue> inputs;
  inputs.push_back(tensor);

  /////////////////////////////////////////////////////////////////////////////
  // Get outputs (all nodes)
  /////////////////////////////////////////////////////////////////////////////
  // All tensors are detached in order to avoid the overhead of gradient calculation and then to
  // save memory/CPU time
  // Tensors are also expected stored in a contiguous block memory

  // Chemical potentials
  auto output_mu_tmp = this->chemical_potentials_nn_[given_phase]->forward(inputs);
  at::Tensor output_mu = output_mu_tmp.toTensor().detach();

  // Mobilities
  at::Tensor output_mob_tensor;
  if (this->has_nn_mobilities_ && this->mobilities_nn_.contains(given_phase)) {
    output_mob_tensor = output_mu;
    if (this->own_mobilities_model_) {
      auto output_mob_tmp = this->mobilities_nn_[given_phase]->forward(inputs);
      output_mob_tensor = output_mob_tmp.toTensor().detach();
    }
    if (!output_mob_tensor.is_contiguous()) {
      output_mob_tensor = output_mob_tensor.contiguous();
    }
  }
  // Energies
  at::Tensor output_nrj_tensor;
  if (this->has_nn_energies_ && this->energies_nn_.contains(given_phase)) {
    output_nrj_tensor = output_mu;
    if (this->own_energies_model_) {
      auto output_nrj_tmp = this->energies_nn_[given_phase]->forward(inputs);
      output_nrj_tensor = output_nrj_tmp.toTensor().detach();
    }
    if (!output_nrj_tensor.is_contiguous()) {
      output_nrj_tensor = output_nrj_tensor.contiguous();
    }
  }

  if (!output_mu.is_contiguous()) {
    output_mu = output_mu.contiguous();
  }

  const int output_mu_size = output_mu.sizes()[1];
  const double* output_data_mu = output_mu.data_ptr<double>();

  int output_mob_size;
  double* output_data_mob;
  if (this->has_nn_mobilities_ && this->mobilities_nn_.contains(given_phase)) {
    output_mob_size = output_mob_tensor.sizes()[1];
    output_data_mob = output_mob_tensor.data_ptr<double>();
  }

  int output_nrj_size;
  double* output_data_nrj;
  if (this->has_nn_energies_ && this->energies_nn_.contains(given_phase)) {
    output_nrj_size = output_nrj_tensor.sizes()[1];
    output_data_nrj = output_nrj_tensor.data_ptr<double>();
  }

  // Fill Calphad containers in order to update outputs of the Calphad problem
  ia = -1;
  for (const auto& i : list_nodes) {
    ia++;
    for (std::size_t iel = 0; iel < chemical_system.size(); ++iel) {
      const std::string& elem = std::get<0>(chemical_system[iel]);

      // Chemical potentials
      const int& mu_index = this->index_chemical_potentials_neural_network_.at(given_phase);
      this->chemical_potentials_[std::make_tuple(i, elem)] =
          output_data_mu[ia * output_mu_size + mu_index + this->idx_elem_[elem]];

      // Mobilities
      if (this->has_nn_mobilities_ && this->mobilities_nn_.contains(given_phase)) {
        const int& mob_index = this->index_mobilities_neural_network_.at(given_phase);
        this->mobilities_[std::make_tuple(i, given_phase, elem)] =
            output_data_mob[ia * output_mob_size + mob_index + this->idx_elem_[elem]];
      }

      // Energies
      if (this->has_nn_energies_ && this->energies_nn_.contains(given_phase)) {
        // ir = index of energy_name
        int ir = 0;
        for (const auto& energy_name : this->input_energies_order_) {
          const auto& key = std::make_tuple(i, given_phase, energy_name);
          const int& energy_index = this->index_energies_neural_network_.at(given_phase);
          this->energies_of_phases_[key] =
              output_data_nrj[ia * output_nrj_size + energy_index + ir];
          ir++;
        }
      }
    }
    // Diffusion chemcial potentials
    for (std::size_t iel = 0; iel < chemical_system.size(); ++iel) {
      const std::string& elem = std::get<0>(chemical_system[iel]);
      if (elem == this->element_removed_from_nn_inputs_) continue;
      this->diffusion_chemical_potentials_[std::make_tuple(i, elem)] =
          this->chemical_potentials_[std::make_tuple(i, elem)] -
          this->chemical_potentials_[std::make_tuple(i, this->element_removed_from_nn_inputs_)];
    }
  }

  // Clear containers
  inputs.clear();
  tp_gf_at_node.clear();
}

/**
 * @brief Check the consistency of the outputs of the Calphad problem.
 *
 * @tparam T
 * @param output_system The outputs of the Calphad problem (primary variables).
 */
template <typename T>
void CalphadInformedNeuralNetwork<T>::check_variables_consistency(
    std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system) {
  for (auto& [output_infos, output_value] : output_system) {
    const std::string& output_type = output_infos.back();

    // Fill output with the relevant values
    switch (calphad_outputs::from(output_type)) {
      case calphad_outputs::cp: {
        MFEM_VERIFY(false,
                    "CalphadInformedNeuralNetwork is only built for mu, mob, h, hm, g, gm.\n");

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
void CalphadInformedNeuralNetwork<T>::finalize() {}

/**
 * @brief Destroy the Binary Melting< T>:: Binary Melting object
 *
 * @tparam T
 */
template <typename T>
CalphadInformedNeuralNetwork<T>::~CalphadInformedNeuralNetwork() {}

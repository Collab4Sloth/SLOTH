/**
 * @file CalphadInformedNeuralNetwork.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Calphad-informed neuronal network used to predict Calphad contributions
 * @version 0.1
 * @date 2025-03-31
 *
 * @copyright Copyright (c) 2025
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

// TODO(cci) : metamodele per phase also for potential
template <typename T>
class CalphadInformedNeuralNetwork : public CalphadBase<T> {
 private:
  std::unique_ptr<torch::jit::script::Module> chemical_potentials_nn_;
  std::unordered_map<std::string, std::unique_ptr<torch::jit::script::Module>> mobilities_nn_;
  std::unique_ptr<torch::jit::script::Module> energies_nn_;
  std::unique_ptr<torch::jit::script::Module> molar_fractions_nn_;

  std::string chemical_potentials_neural_network_;
  std::unordered_map<std::string, std::string> mobilities_neural_network_;
  std::string energies_neural_network_;
  std::string molar_fractions_neural_network_;

  int index_chemical_potentials_neural_network_;
  std::unordered_map<std::string, int> index_mobilities_neural_network_;
  int index_energies_neural_network_;
  int index_molar_fractions_neural_network_;

  std::map<std::string, unsigned int> idx_elem_;
  std::vector<std::string> input_composition_order_;
  double input_composition_factor_;

  bool model_built_with_pressure_{false};
  bool own_mobilities_model_{false};

  std::unique_ptr<CalphadUtils<T>> CU_;
  void compute(const std::set<int>& list_nodes, const std::vector<T>& tp_gf,
               const std::vector<std::tuple<std::string, std::string>>& chemical_system);

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

  // Chemical potentials
  if (this->params_.has_parameter("ChemicalPotentialsNeuralNetwork")) {
    this->chemical_potentials_neural_network_ =
        this->params_.template get_param_value<std::string>("ChemicalPotentialsNeuralNetwork");
    check_file_existence(this->chemical_potentials_neural_network_);
    this->index_chemical_potentials_neural_network_ =
        this->params_.template get_param_value_or_default<int>(
            "ChemicalPotentialsNeuralNetworkIndex", 0);
  }

  // Mobilities
  if (this->params_.has_parameter("MobilitiesNeuralNetwork")) {
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

  // Check if a specific model is used for mobilities. Otherwise the model for chemical potential is
  // used
  if (this->params_.has_parameter("OwnMobilityModel")) {
    this->own_mobilities_model_ = this->params_.template get_param_value<bool>("OwnMobilityModel");
  }
}

////////////////////////////////
////////////////////////////////
/**
 * @brief Construct a new CalphadInformedNeuralNetwork::CalphadInformedNeuralNetwork object
 *
 * @param params
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
 * @param params
 * @param is_KKS
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
      this->chemical_potentials_nn_ = std::make_unique<torch::jit::script::Module>(
          torch::jit::load(this->chemical_potentials_neural_network_));
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
void CalphadInformedNeuralNetwork<T>::execute(
    const int dt, const std::set<int>& list_nodes, const std::vector<T>& tp_gf,
    const std::vector<std::tuple<std::string, std::string>>& chemical_system,
    std::optional<std::vector<std::tuple<std::string, std::string, double>>> status_phase) {
  Catch_Time_Section("CalphadInformedNeuralNetwork<T>::execute");
  // Clear containers and recalculation of the numbers of nodes
  // const size_t nb_nodes = this->CU_->get_size(tp_gf[0]);
  // this->clear_containers();

  // Thermodynamic Calculations
  this->compute(list_nodes, tp_gf, chemical_system);

  // Use containers to update output_system
  // this->update_outputs(nb_nodes, output_system);
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
    const std::vector<std::tuple<std::string, std::string>>& chemical_system) {
  Catch_Time_Section("CalphadInformedNeuralNetwork<T>::compute");
  const size_t nb_nodes = list_nodes.size();
  std::vector<double> tp_gf_at_node(tp_gf.size());

  // TODO(cci) add pressure as default inputs
  int nb_input = chemical_system.size() + 1;
  if (this->model_built_with_pressure_) {
    nb_input++;
  }
  auto options = torch::TensorOptions().dtype(torch::kDouble);
  auto tensor =
      torch::zeros({static_cast<int64_t>(nb_nodes), static_cast<int64_t>(nb_input)}, options);

  for (const auto& i : list_nodes) {
    // Populate tp_gf_at_node for the current node
    std::transform(tp_gf.begin(), tp_gf.end(), tp_gf_at_node.begin(),
                   [&i](const T& vec) { return vec[i]; });

    const double temperature = tp_gf_at_node[0];
    const double pressure = tp_gf_at_node[1];
    tensor[i][0] = temperature;

    // id_elem starts to 1 if pressure is not an input of the model. Otherwise it starts to 2.
    int id_elem = 1;
    if (this->model_built_with_pressure_) {
      tensor[i][1] = pressure;
      id_elem = 2;
    }
    //
    for (std::size_t iel = 0; iel < chemical_system.size(); ++iel) {
      const std::string& elem = std::get<0>(chemical_system[iel]);
      tensor[i][this->idx_elem_[elem] + id_elem] =
          tp_gf_at_node[iel + 2] * input_composition_factor_;
    }
  }

  std::vector<torch::jit::IValue> inputs;
  inputs.push_back(tensor);
  // Get CALPHAD contributions from Neural Networks

  std::unordered_map<std::string, at::Tensor> output_mob;
  std::unordered_map<std::string, int> output_mob_size;
  std::unordered_map<std::string, const double*> output_data_mob;
  auto output_mu_tmp = this->chemical_potentials_nn_->forward(inputs);

  at::Tensor output_mu = output_mu_tmp.toTensor().detach();

  for (const auto& [phase, mob_model] : this->mobilities_nn_) {
    at::Tensor output_mob_tensor = output_mu;
    if (this->own_mobilities_model_) {
      output_mob_tensor = mob_model->forward(inputs).toTensor().detach();
    }
    if (!output_mob_tensor.is_contiguous()) {
      output_mob_tensor = output_mob_tensor.contiguous();
    }
    output_mob[phase] = output_mob_tensor;
    output_mob_size[phase] = output_mob_tensor.sizes()[1];
    output_data_mob[phase] = output_mob_tensor.data_ptr<double>();
  }

  // auto output_mu_tmp = this->chemical_potentials_nn_->forward(inputs);

  // at::Tensor output_mu = output_mu_tmp.toTensor().detach();
  if (!output_mu.is_contiguous()) {
    output_mu = output_mu.contiguous();
  }
  const int output_mu_size = output_mu.sizes()[1];
  const double* output_data_mu = output_mu.data_ptr<double>();
  for (const auto& i : list_nodes) {
    for (std::size_t iel = 0; iel < chemical_system.size(); ++iel) {
      const std::string& elem = std::get<0>(chemical_system[iel]);
      this->chemical_potentials_[std::make_tuple(i, elem)] =
          output_data_mu[i * output_mu_size + this->index_chemical_potentials_neural_network_ +
                         this->idx_elem_[elem]];
      // Mobility of element iel for each phase
      for (auto& [phase, output] : output_data_mob) {
        const int& mindex = this->index_mobilities_neural_network_.at(phase);
        const int& msize = output_mob_size.at(phase);
        this->mobilities_[std::make_tuple(i, phase, elem)] =
            output[i * msize + mindex + this->idx_elem_[elem]];
      }
    }
  }

  // Clear containers
  inputs.clear();
  tp_gf_at_node.clear();
  output_mob.clear();
  output_mob_size.clear();
  output_data_mob.clear();
}

/**
 * @brief Check the consistency of outputs required for the current Calphad problem
 *
 * @tparam T
 * @param output_system
 */
template <typename T>
void CalphadInformedNeuralNetwork<T>::check_variables_consistency(
    std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system) {
  for (auto& [output_infos, output_value] : output_system) {
    const std::string& output_type = output_infos.back();

    // Fill output with the relevant values
    switch (calphad_outputs::from(output_type)) {
      case calphad_outputs::g:
      case calphad_outputs::gm:
      case calphad_outputs::h:
      case calphad_outputs::hm:
      case calphad_outputs::cp: {
        MFEM_VERIFY(false, "CalphadInformedNeuralNetwork is only built for mu, mob.\n");

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

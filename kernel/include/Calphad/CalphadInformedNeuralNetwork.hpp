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
#pragma once

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

#include "Calphad/CalphadInformedNeuralNetwork.tpp"

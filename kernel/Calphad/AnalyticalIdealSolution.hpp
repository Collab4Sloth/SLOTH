/**
 * @file AnalyticalIdealSolution.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Analytical thermodynamic description for an ideal solution
 * @version 0.1
 * @date 2025-01-07
 *
 * Copyright CEA (c) 2025
 *
 */
#include <algorithm>
#include <functional>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "Calphad/CalphadBase.hpp"
#include "Calphad/CalphadUtils.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"

#pragma once

template <typename T>
class AnalyticalIdealSolution : public CalphadBase<T> {
 private:
  std::unique_ptr<CalphadUtils<T>> CU_;
  void compute(const std::set<int>& list_nodes, const std::vector<T>& tp_gf,
               const std::vector<std::tuple<std::string, std::string>>& chemical_system);

  void check_variables_consistency(
      std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system);

 public:
  constexpr explicit AnalyticalIdealSolution(const Parameters& params);
  constexpr AnalyticalIdealSolution(const Parameters& params, bool is_KKS);

  void initialize(
      const std::vector<std::tuple<std::string, std::string>>& sorted_chemical_system) override;

  void execute(const int dt, const std::set<int>& list_nodes, const std::vector<T>& aux_gf,
               const std::vector<std::tuple<std::string, std::string>>& chemical_system,
               std::optional<std::vector<std::tuple<std::string, std::string>>> status_phase =
                   std::nullopt) override;

  void finalize() override;

  ////////////////////////////////

  void get_parameters() override;

  ////////////////////////////////

  ~AnalyticalIdealSolution();
};

////////////////////////////////
////////////////////////////////

/**
 * @brief Get the parameters associated with the AnalyticalIdealSolution object
 *
 * @tparam T
 */
template <typename T>
void AnalyticalIdealSolution<T>::get_parameters() {
  CalphadBase<T>::get_parameters();
  this->description_ = this->params_.template get_param_value_or_default<std::string>(
      "description", "Analytical thermodynamic description for an ideal solution ");
}

////////////////////////////////
////////////////////////////////
/**
 * @brief Construct a new AnalyticalIdealSolution::AnalyticalIdealSolution object
 *
 * @param params
 */
template <typename T>
constexpr AnalyticalIdealSolution<T>::AnalyticalIdealSolution(const Parameters& params)
    : CalphadBase<T>(params, false) {
  this->CU_ = std::make_unique<CalphadUtils<T>>();
  this->get_parameters();
}

/**
 * @brief Construct a new AnalyticalIdealSolution::AnalyticalIdealSolution object
 *
 * @param params
 * @param is_KKS
 */
template <typename T>
constexpr AnalyticalIdealSolution<T>::AnalyticalIdealSolution(const Parameters& params, bool is_KKS)
    : CalphadBase<T>(params, is_KKS) {
  this->CU_ = std::make_unique<CalphadUtils<T>>();
  this->get_parameters();
}

/**
 * @brief Initialization of the thermodynamic calculation
 *
 * @tparam T
 */
template <typename T>
void AnalyticalIdealSolution<T>::initialize(
    const std::vector<std::tuple<std::string, std::string>>& sorted_chemical_system) {}

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
void AnalyticalIdealSolution<T>::execute(
    const int dt, const std::set<int>& list_nodes, const std::vector<T>& tp_gf,
    const std::vector<std::tuple<std::string, std::string>>& chemical_system,
    std::optional<std::vector<std::tuple<std::string, std::string>>> status_phase) {
  // Clear containers and recalculation of the numbers of nodes
  // const size_t nb_nodes = this->CU_->get_size(tp_gf[0]);
  // this->clear_containers();
  // if (dt == 1) {
  //   this->check_variables_consistency(output_system);
  // }
  // Thermodynamic Calculations
  this->compute(list_nodes, tp_gf, chemical_system);

  // Use containers to update output_system
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
void AnalyticalIdealSolution<T>::compute(
    const std::set<int>& list_nodes, const std::vector<T>& tp_gf,
    const std::vector<std::tuple<std::string, std::string>>& chemical_system) {
  // Let us assume an ideal mixing solution
  const std::string& phase = "SOLUTION";
  const std::vector<std::string> energy_names = {"G"};
  std::vector<double> tp_gf_at_node(tp_gf.size());

  const auto temperature_sort_method =
      this->params_.template get_param_value_or_default<std::string>("temperature_sort_method",
                                                                     "No");
  const auto pressure_sort_method =
      this->params_.template get_param_value_or_default<std::string>("pressure_sort_method", "No");

  std::vector<int> sorted_n_t_p =
      this->CU_->sort_nodes(tp_gf[0], tp_gf[1], temperature_sort_method, pressure_sort_method);

  // Remove elements from list_nodes in sorted_n_t_p
  // Usefull if the targeted list is lower than the complete list (eg. KKS problems)
  sorted_n_t_p.erase(
      std::remove_if(sorted_n_t_p.begin(), sorted_n_t_p.end(),
                     [&list_nodes](int node) { return list_nodes.find(node) != list_nodes.end(); }),
      sorted_n_t_p.end());

  // Process CALPHAD calculations for each node
  for (const auto& id : sorted_n_t_p) {
    // Populate tp_gf_at_node for the current node
    std::transform(tp_gf.begin(), tp_gf.end(), tp_gf_at_node.begin(),
                   [&id](const T& vec) { return vec[id]; });
    const auto Temp = tp_gf_at_node[0];
    const auto& epsilon = 1.e-10;
    const auto x = std::max(epsilon, std::min(1. - epsilon, tp_gf_at_node[2]));
    const auto& elem = std::get<0>(chemical_system[0]);

    // Energy
    for (const auto& energy_name : energy_names) {
      const auto& key = std::make_tuple(id, phase, energy_name);
      this->energies_of_phases_[key] =
          Physical::R * Temp * (x * std::log(x) + (1. - x) * std::log(1. - x));
    }

    // Chemical potential
    this->chemical_potentials_[std::make_tuple(id, elem)] =
        Physical::R * Temp * (std::log(x) - std::log(1. - x));

    // Molar fraction
    const auto& key = std::make_tuple(id, phase, elem);
    this->elem_mole_fraction_by_phase_[key] = x;
  }
}
/**
 * @brief Check the consistency of outputs required for the current Calphad problem
 *
 * @tparam T
 * @param output_system
 */
template <typename T>
void AnalyticalIdealSolution<T>::check_variables_consistency(
    std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system) {
  for (auto& [output_infos, output_value] : output_system) {
    const std::string& output_type = output_infos.back();

    // Fill output with the relevant values
    switch (calphad_outputs::from(output_type)) {
      case calphad_outputs::gm:
      case calphad_outputs::h:
      case calphad_outputs::hm: {
        MFEM_VERIFY(false, "AnalyticalIdealSolution is only built for mu, x and g.\n");

        SlothInfo::debug("Output not available for this Calphad problem: ", output_type);
        break;
      }
      // (Optional) handle allowed cases explicitly
      case calphad_outputs::mu:
      case calphad_outputs::x:
      case calphad_outputs::g:
        // process normally
        break;

      default:
        MFEM_VERIFY(false, "Unhandled output type.\n");
        break;
    }
  }
}

/**
 * @brief Finalization actions (free memory)
 *
 * @tparam T
 */
template <typename T>
void AnalyticalIdealSolution<T>::finalize() {}

/**
 * @brief Destroy the Binary Melting< T>:: Binary Melting object
 *
 * @tparam T
 */
template <typename T>
AnalyticalIdealSolution<T>::~AnalyticalIdealSolution() {}

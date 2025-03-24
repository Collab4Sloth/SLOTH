/**
 * @file BinaryTabulation.hpp
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
#include <memory>
#include <string>
#include <tuple>
#include <vector>
#include <utility>

#include "Calphad/CalphadBase.hpp"
#include "Calphad/CalphadUtils.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"

#pragma once

template <typename T>
class BinaryTabulation : CalphadBase<T> {
 private:
  std::vector<double> tab_x;
  std::vector<double> tab_g;
  std::vector<double> tab_m;
  std::unique_ptr<CalphadUtils<T>> CU_;
  void compute(
      const size_t nb_nodes, const std::vector<T>& tp_gf,
      const std::vector<std::tuple<std::string, std::string>>& chemical_system,
      std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system);

  void check_variables_consistency(
      std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system);

 public:
  explicit BinaryTabulation(const Parameters& params);

  void initialize() override;

  void execute(const int dt, const std::vector<T>& aux_gf,
               const std::vector<std::tuple<std::string, std::string>>& chemical_system,
               std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>&
                   output_system) override;

  void finalize() override;

  ////////////////////////////////

  void get_parameters() override;

  ////////////////////////////////

  ~BinaryTabulation();
};

////////////////////////////////
////////////////////////////////

/**
 * @brief Get the parameters associated with the BinaryTabulation object
 *
 * @tparam T
 */
template <typename T>
void BinaryTabulation<T>::get_parameters() {
  this->description_ = this->params_.template get_param_value_or_default<std::string>(
      "description",
      "Analytical thermodynamic description for an ideal solution using tabulated data ");
}

////////////////////////////////
////////////////////////////////
/**
 * @brief Construct a new BinaryTabulation::BinaryTabulation object
 *
 * @param params
 */
template <typename T>
BinaryTabulation<T>::BinaryTabulation(const Parameters& params) : CalphadBase<T>(params) {
  this->CU_ = std::make_unique<CalphadUtils<T>>();
  this->get_parameters();
}

/**
 * @brief Initialization of the thermodynamic calculation
 *
 * @tparam T
 */
template <typename T>
void BinaryTabulation<T>::initialize() {
  std::string filename = this->params_.template get_param_value_or_default<std::string>(
      "data filename", "no input file");
  std::ifstream file(filename);
  MFEM_VERIFY(file, "Error: file " << filename << " does not exist");

  std::string line;

  while (std::getline(file, line)) {
    std::stringstream ss(line);
    std::string value;
    std::vector<double> row;

    while (std::getline(ss, value, ',')) {
      row.emplace_back(std::stod(value));
    }

    if (!row.empty()) {
      tab_x.emplace_back(std::move(row[0]));
      tab_g.emplace_back(std::move(row[1]));
      tab_m.emplace_back(std::move(row[2]));
    }
  }
  file.close();
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
void BinaryTabulation<T>::execute(
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
void BinaryTabulation<T>::compute(
    const size_t nb_nodes, const std::vector<T>& tp_gf,
    const std::vector<std::tuple<std::string, std::string>>& chemical_system,
    std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system) {
  // Let us assume an ideal mixing solution
  const std::string& phase = "SOLUTION";
  const std::vector<std::string> energy_names = {"G"};
  std::vector<double> tp_gf_at_node(tp_gf.size());

  std::vector<int> sorted_n_t_p = this->CU_->sort_nodes(tp_gf[0], tp_gf[1], "No", "No");

  // Process CALPHAD calculations for each node
  for (const auto& id : sorted_n_t_p) {
    // Populate tp_gf_at_node for the current node
    std::transform(tp_gf.begin(), tp_gf.end(), tp_gf_at_node.begin(),
                   [&id](const T& vec) { return vec[id]; });
    const auto Temp = tp_gf_at_node[0];
    const auto& epsilon = 1.e-10;
    int i = 0;
    const auto& elem = std::get<0>(chemical_system[0]);
    const auto x = tp_gf_at_node[2];

    auto lower = std::lower_bound(tab_x.begin(), tab_x.end(), x);
    if (lower != tab_x.begin() && (lower + 1) != tab_x.end()) {
      i = lower - tab_x.begin() - 1;
    } else if ((lower + 1) == tab_x.end()) {
      i = tab_x.size() - 1;
    }
    // Energy
    for (const auto& energy_name : energy_names) {
      const auto& key = std::make_tuple(id, phase, energy_name);

      this->energies_of_phases_[key] =
          tab_g[i] + (tab_g[i + 1] - tab_g[i]) * (x - tab_x[i]) / (tab_x[i + 1] - tab_x[i]);
    }

    // Chemical potential
    this->chemical_potentials_[std::make_tuple(id, elem)] =
        tab_m[i] + (tab_m[i + 1] - tab_m[i]) * (x - tab_x[i]) / (tab_x[i + 1] - tab_x[i]);

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
void BinaryTabulation<T>::check_variables_consistency(
    std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system) {
  for (auto& [output_infos, output_value] : output_system) {
    const std::string& output_type = output_infos.back();

    // Fill output with the relevant values
    switch (calphad_outputs::from(output_type)) {
      case calphad_outputs::gm:
      case calphad_outputs::h:
      case calphad_outputs::hm: {
        MFEM_VERIFY(false, "BinaryTabulation is only built for mu, x and g.\n");

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
void BinaryTabulation<T>::finalize() {}

/**
 * @brief Destroy the Binary Melting< T>:: Binary Melting object
 *
 * @tparam T
 */
template <typename T>
BinaryTabulation<T>::~BinaryTabulation() {}

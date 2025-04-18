
/**
 * @file CalphadBase.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Base class for Calphad objets
 * @version 0.1
 * @date 2025-01-07
 *
 * Copyright CEA (c) 2025
 *
 */
#include <algorithm>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

#include "Calphad/CalphadUtils.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Profiling/Profiling.hpp"

#pragma once

template <typename T>
class CalphadBase {
 private:
  bool is_KKS_;
  std::string KKS_secondary_phase_;
  double KKS_temperature_increment_;
  double KKS_composition_increment_;
  double KKS_threshold_;
  std::string KKS_temperature_scheme_;

  void get_KKS_parameters();
  void KKS_execute(const int dt, const std::vector<T> &tp_gf,
                   const std::tuple<std::string, T> &phasefields_gf,
                   const std::vector<std::tuple<std::string, std::string>> &chemicalsystem);

 protected:
  std::unique_ptr<CalphadUtils<T>> CU_;
  std::string description_{"UNKNOWN CALPHAD"};

  const Parameters &params_;
  MapStringDouble unsuspended_phases_;

  // Containers used to store the results of the equilibrium calculations
  // Chemical potentials for each element. Nodal values

  // Chemical potential. Nodal values
  // key: [node, elem]
  std::map<std::tuple<int, std::string>, double> chemical_potentials_;
  std::map<std::tuple<int, std::string, std::string>, double> chemical_potentials_left_T_;
  std::map<std::tuple<int, std::string, std::string>, double> chemical_potentials_left_x_;
  std::map<std::tuple<int, std::string, std::string>, double> chemical_potentials_right_T_;
  std::map<std::tuple<int, std::string, std::string>, double> chemical_potentials_right_x_;
  // Element molar fraction for each phase .Nodal values
  // Element mole fraction by phase. Nodal values
  // key: [node, phase, elem]
  std::map<std::tuple<int, std::string, std::string>, double> elem_mole_fraction_by_phase_;
  // Site fraction of a given constituant, in a given sublattice for each phase. Nodal values
  // Site fraction. Nodal values
  // key: [node, phase, cons, sub]
  std::map<std::tuple<int, std::string, std::string, int>, double> site_fraction_;
  // Energy for each phase. Nodal values
  // key: [node, phase, symbol]
  // Energy of each phase. Nodal values
  // key: [node, phase, energy_type]
  std::map<std::tuple<int, std::string, std::string>, double> energies_of_phases_;
  // Driving force for each phase. Nodal values
  // key: [node, phase]
  std::map<std::tuple<int, std::string>, double> driving_forces_;
  // Mobility for each element in a given phase. Nodal values
  // key: [node, phase, elem]
  std::map<std::tuple<int, std::string, std::string>, double> mobilities_;
  //  Heat capacity. Nodal values
  // key: [node]
  std::map<int, double> heat_capacity_;

  void clear_containers();

  void update_outputs(
      const size_t nb_nodes,
      std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>> &output_system);

  virtual void execute(const int dt, const std::set<int> &list_nodes, const std::vector<T> &tp_gf,
                       const std::vector<std::tuple<std::string, std::string>> &chemicalsystem,
                       std::optional<std::string> phase = std::nullopt) = 0;

 public:
  constexpr explicit CalphadBase(const Parameters &params);
  constexpr CalphadBase(const Parameters &params, bool is_KKS);

  virtual void initialize(
      const std::vector<std::tuple<std::string, std::string>> &sorted_chemical_system) = 0;

  void global_execute(
      const int dt, const std::vector<T> &tp_gf,
      const std::vector<std::tuple<std::string, std::string>> &chemicalsystem,
      std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>> &output_system,
      std::optional<const std::tuple<std::string, T>> phase_fields = std::nullopt);

  virtual void finalize() = 0;

  ////////////////////////////////

  ////////////////////////////////

  virtual void get_parameters();

  ////////////////////////////////

  virtual ~CalphadBase();

  std::string get_description() { return this->description_; }
};

/**
 * @brief Construct a new Calphad Base< T>:: Calphad Base object
 *
 * @tparam T
 */
template <typename T>
constexpr CalphadBase<T>::CalphadBase(const Parameters &params) : CalphadBase(params, false) {
  this->CU_ = std::make_unique<CalphadUtils<T>>();
  this->get_parameters();
}

/**
 * @brief Construct a new Calphad Base< T>:: Calphad Base object
 *
 * @tparam T
 * @param params
 * @param is_KKS
 */
template <typename T>
constexpr CalphadBase<T>::CalphadBase(const Parameters &params, bool is_KKS)
    : params_(params), is_KKS_(is_KKS) {
  this->CU_ = std::make_unique<CalphadUtils<T>>();

  this->get_parameters();
}

/**
 * @brief Get Calphad parameters
 *
 * @tparam T
 */
template <typename T>
void CalphadBase<T>::get_parameters() {
  if (this->is_KKS_) {
    this->get_KKS_parameters();
  }
}
/**
 * @brief Get all parameters required by KKS problem
 *
 * @tparam T
 */
template <typename T>
void CalphadBase<T>::get_KKS_parameters() {
  this->KKS_secondary_phase_ =
      this->params_.template get_param_value<std::string>("KKS_secondary_phase");

  this->KKS_temperature_increment_ =
      this->params_.template get_param_value<double>("KKS_temperature_increment");
  this->KKS_composition_increment_ =
      this->params_.template get_param_value<double>("KKS_composition_increment");
  // this->KKS_temperature_scheme_ =
  //     this->params_.template get_param_value<std::string>("KKS_temperature_scheme");
  this->KKS_threshold_ =
      this->params_.template get_param_value_or_default<double>("KKS_threshold", 1.e-2);
}

/**
 * @brief Compute Calphad contributions needed by KKS resolution
 *
 * @tparam T
 * @param dt
 * @param tp_gf
 * @param chemicalsystem
 */
template <typename T>
void CalphadBase<T>::KKS_execute(
    const int dt, const std::vector<T> &tp_gf, const std::tuple<std::string, T> &phasefields_gf,
    const std::vector<std::tuple<std::string, std::string>> &chemicalsystem) {
  // Creation initial list of nodes
  const size_t nb_nodes = this->CU_->get_size(tp_gf[0]);
  std::set<int> list_nodes;
  for (int i = 0; i < nb_nodes; ++i) {
    list_nodes.insert(i);
  }
  ////////////////////////////////////////////////////////
  // List of nodes by phases and within interface
  ////////////////////////////////////////////////////////
  const auto &[phase, phi_gf] = phasefields_gf;

  std::set<int> indices_ph_1;
  std::set<int> indices_ph_2;
  std::set<int> indices_inter;
  for (int index = 0; index < nb_nodes; ++index) {
    if (phi_gf[index] > 1 - this->KKS_threshold_) {
      indices_ph_1.insert(index);
    } else if (phi_gf[index] < this->KKS_threshold_) {
      indices_ph_2.insert(index);
    } else {
      indices_inter.insert(index);
    }
  }
  ////////////////////////////////////////////////////////
  /// Calculation of all thermodynamic contribution
  ////////////////////////////////////////////////////////
  // Primary phase
  this->execute(dt, indices_ph_1, tp_gf, chemicalsystem, phase);
  // Secondary phase
  this->execute(dt, indices_ph_2, tp_gf, chemicalsystem, this->KKS_secondary_phase_);

  // Interface calculations
  auto calculate_interface = [&](const std::vector<T> &delta_tp_gf, const std::string &given_phase,
                                 std::map<std::tuple<int, std::string, std::string>, double>
                                     &chemical_potential_interface) {
    this->execute(dt, indices_inter, delta_tp_gf, chemicalsystem, given_phase);
    for (const auto &in : indices_inter) {
      for (const auto &elem : chemicalsystem) {
        const auto &[elem1, unit] = elem;
        const double mu = this->chemical_potentials_.at(std::make_tuple(in, elem1));
        chemical_potential_interface.emplace(std::make_tuple(in, elem1, given_phase), mu);
      }
    }
  };

  // T-dT
  std::vector<T> delta_tp_gf = tp_gf;
  delta_tp_gf[0] += this->KKS_temperature_increment_;
  calculate_interface(delta_tp_gf, phase, this->chemical_potentials_left_T_);
  calculate_interface(delta_tp_gf, this->KKS_secondary_phase_, this->chemical_potentials_left_T_);

  // T+dT
  delta_tp_gf = tp_gf;
  delta_tp_gf[0] -= this->KKS_temperature_increment_;
  calculate_interface(delta_tp_gf, phase, this->chemical_potentials_right_T_);
  calculate_interface(delta_tp_gf, this->KKS_secondary_phase_, this->chemical_potentials_right_T_);

  // Loop over chemicalsystem
  for (std::size_t i = 0; i < chemicalsystem.size(); ++i) {
    // x+dx
    delta_tp_gf = tp_gf;
    delta_tp_gf[i + 2] += this->KKS_composition_increment_;
    calculate_interface(delta_tp_gf, phase, this->chemical_potentials_right_x_);
    calculate_interface(delta_tp_gf, this->KKS_secondary_phase_,
                        this->chemical_potentials_right_x_);

    // x-dx
    delta_tp_gf = tp_gf;
    delta_tp_gf[i + 2] -= this->KKS_composition_increment_;
    calculate_interface(delta_tp_gf, phase, this->chemical_potentials_left_x_);
    calculate_interface(delta_tp_gf, this->KKS_secondary_phase_, this->chemical_potentials_left_x_);
  }
  ////////////////////////////////////////////////////////
  /// Solve linear systems
  ////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////
  /// Recover thermodynamic contribution
  ////////////////////////////////////////////////////////
}

/**
 * @brief High-level method for managing equilibrium calculations
 *
 * @tparam T
 * @param dt
 * @param tp_gf
 * @param chemicalsystem
 * @param output_system
 */
template <typename T>
void CalphadBase<T>::global_execute(
    const int dt, const std::vector<T> &tp_gf,
    const std::vector<std::tuple<std::string, std::string>> &chemicalsystem,
    std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>> &output_system,
    std::optional<const std::tuple<std::string, T>> phase_field_gf) {
  const size_t nb_nodes = this->CU_->get_size(tp_gf[0]);
  // Reinitialize containers
  this->clear_containers();
  // Execute
  if (!this->is_KKS_) {
    // Creation list of nodes
    std::set<int> list_nodes;
    for (int i = 0; i < nb_nodes; ++i) {
      list_nodes.insert(i);
    }
    this->execute(dt, list_nodes, tp_gf, chemicalsystem);
    list_nodes.clear();
  } else {
    MFEM_VERIFY(phase_field_gf.has_value(),
                "Error: phase_fields_gf is required for KKS execution.");
    this->KKS_execute(dt, tp_gf, *phase_field_gf, chemicalsystem);
  }
  // Use containers to update output_system
  this->update_outputs(nb_nodes, output_system);
}

/**
 * @brief Clear containers used to store the results of equilibrium calculations
 *
 * @tparam T
 */
template <typename T>
void CalphadBase<T>::clear_containers() {
  // Clear before filling with new results
  this->chemical_potentials_.clear();
  this->elem_mole_fraction_by_phase_.clear();
  this->site_fraction_.clear();
  this->energies_of_phases_.clear();
  this->driving_forces_.clear();
  this->heat_capacity_.clear();
  this->mobilities_.clear();
}

/**
 * @brief Update the outputs on the basis of results stored in containers
 *
 * @tparam T
 * @param nb_nodes
 * @param output_system
 */
template <typename T>
void CalphadBase<T>::update_outputs(
    const size_t nb_nodes,
    std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>> &output_system) {
  Catch_Time_Section("CalphadBase<T>::update_outputs");

  ////////////////////
  // Update outputs //
  ////////////////////
  T output(nb_nodes);

  for (auto &[output_infos, output_value] : output_system) {
    const std::string &output_type = output_infos.back();

    // Fill output with the relevant values
    switch (calphad_outputs::from(output_type)) {
      case calphad_outputs::mu: {
        const std::string &output_elem = output_infos[1];
        for (std::size_t i = 0; i < nb_nodes; ++i) {
          output[i] = this->chemical_potentials_[std::make_tuple(i, output_elem)];
        }
        break;
      }
      case calphad_outputs::x: {
        const std::string &output_elem = output_infos[1];
        const std::string &output_phase = output_infos[2];
        for (std::size_t i = 0; i < nb_nodes; ++i) {
          output[i] =
              this->elem_mole_fraction_by_phase_[std::make_tuple(i, output_phase, output_elem)];
        }
        break;
      }
      case calphad_outputs::y: {
        const std::string &output_cons = output_infos[1];
        const int &output_sub = std::stoi(output_infos[2]);
        const std::string &output_phase = output_infos[3];
        for (std::size_t i = 0; i < nb_nodes; ++i) {
          output[i] =
              this->site_fraction_[std::make_tuple(i, output_phase, output_cons, output_sub)];
        }
        break;
      }
      case calphad_outputs::g: {
        const std::string &output_phase = output_infos[1];
        for (std::size_t i = 0; i < nb_nodes; ++i) {
          output[i] = this->energies_of_phases_[std::make_tuple(i, output_phase, "G")];
        }
        break;
      }
      case calphad_outputs::gm: {
        const std::string &output_phase = output_infos[1];
        for (std::size_t i = 0; i < nb_nodes; ++i) {
          output[i] = this->energies_of_phases_[std::make_tuple(i, output_phase, "GM")];
        }
        break;
      }
      case calphad_outputs::h: {
        const std::string &output_phase = output_infos[1];
        for (std::size_t i = 0; i < nb_nodes; ++i) {
          output[i] = this->energies_of_phases_[std::make_tuple(i, output_phase, "H")];
        }
        break;
      }
      case calphad_outputs::hm: {
        const std::string &output_phase = output_infos[1];
        for (std::size_t i = 0; i < nb_nodes; ++i) {
          output[i] = this->energies_of_phases_[std::make_tuple(i, output_phase, "HM")];
        }
        break;
      }
      case calphad_outputs::dgm: {
        const std::string &output_phase = output_infos[1];
        for (std::size_t i = 0; i < nb_nodes; ++i) {
          output[i] = this->driving_forces_[std::make_tuple(i, output_phase)];
        }
        break;
      }
      case calphad_outputs::cp: {
        for (std::size_t i = 0; i < nb_nodes; ++i) {
          output[i] = this->heat_capacity_[i];
        }
        break;
      }
      case calphad_outputs::mob: {
        const std::string &output_phase = output_infos[1];
        const std::string &output_elem = output_infos[2];
        for (std::size_t i = 0; i < nb_nodes; ++i) {
          output[i] = this->mobilities_[std::make_tuple(i, output_phase, output_elem)];
        }
        break;
      }
    }
    // Update the referenced output vector
    output_value.get() = output;
  }
}

/**
 * @brief Destroy the Calphad Base< T>:: Calphad Base object
 *
 * @tparam T
 */
template <typename T>
CalphadBase<T>::~CalphadBase() {}

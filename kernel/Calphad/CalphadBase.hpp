
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
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "Calphad/CalphadUtils.hpp"
#include "Coefficients/PhaseFieldPotentials.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Profiling/Profiling.hpp"

#pragma once
// Previous declaration of the KKS class
template <typename T>
class KKS;

template <typename T>
class CalphadBase {
 private:
  bool is_KKS_;

 protected:
  std::shared_ptr<KKS<T>> KKS_;

  std::shared_ptr<CalphadUtils<T>> CU_;

  void clear_containers();

  void update_outputs(
      const int dt, const size_t nb_nodes,
      std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>> &output_system,
      const std::vector<std::tuple<std::vector<std::string>, mfem::Vector>>
          &previous_output_system);
  // Mobility for each element in a given phase. Nodal values. key: [node, phase, elem]
  std::map<std::tuple<int, std::string, std::string>, double> mobilities_;

  // Site fraction of a given constituant, in a given sublattice for each phase. Nodal values
  // Site fraction. Nodal values. key: [node, phase, cons, sub]
  std::map<std::tuple<int, std::string, std::string, int>, double> site_fraction_;

  //  Heat capacity. Nodal values. key: [node]
  std::map<int, double> heat_capacity_;

 public:
  const Parameters &params_;
  std::string element_removed_from_ic_;

  // Containers used to store the results of the equilibrium calculations
  std::map<int, int> error_equilibrium_;

  // Chemical potential. Nodal values. key: [node, elem]
  std::map<std::tuple<int, std::string>, double> chemical_potentials_;
  std::map<std::tuple<int, std::string>, double> diffusion_chemical_potentials_;

  // Mole fraction of phase. Nodal values. key: [node, phase]
  std::map<std::tuple<int, std::string>, double> mole_fraction_of_phase_;

  // Element mole fraction by phase. Nodal values. key: [node, phase, elem]
  std::map<std::tuple<int, std::string, std::string>, double> elem_mole_fraction_by_phase_;

  // Energy for each phase. Nodal values. key: [node, phase, symbol]
  std::map<std::tuple<int, std::string, std::string>, double> energies_of_phases_;

  // Driving force for each phase. Nodal values. key: [node, phase]
  std::map<std::tuple<int, std::string>, double> driving_forces_;
  std::map<std::tuple<int, std::string>, double> nucleus_;

  explicit CalphadBase(const Parameters &params);
  CalphadBase(const Parameters &params, bool is_KKS);

  virtual void initialize(
      const std::vector<std::tuple<std::string, std::string>> &sorted_chemical_system) = 0;

  void global_execute(
      const int dt, const double time_step, const std::vector<T> &tp_gf,
      const std::vector<std::tuple<std::string, std::string>> &chemicalsystem,
      std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>> &output_system,
      const std::vector<std::tuple<std::vector<std::string>, mfem::Vector>> &previous_output_system,
      std::optional<const std::tuple<std::string, T, T>> phase_fields = std::nullopt,
      std::optional<const std::vector<T>> tp_gf_old = std::nullopt,
      std::optional<const std::vector<std::tuple<std::string, std::string, T, T>>> x_gf =
          std::nullopt,
      std::optional<const std::vector<std::tuple<std::string, std::vector<double>>>> coordinates =
          std::nullopt);

  virtual void execute(const int dt, const std::set<int> &list_nodes, const std::vector<T> &tp_gf,
                       const std::vector<std::tuple<std::string, std::string>> &chemicalsystem,
                       std::optional<std::vector<std::tuple<std::string, std::string, double>>>
                           status_phase = std::nullopt) = 0;

  virtual void finalize() = 0;

  ////////////////////////////////

  ////////////////////////////////

  virtual void get_parameters();

  ////////////////////////////////

  virtual ~CalphadBase();
};

/**
 * @brief Construct a new CalphadBase object.
 *
 * @tparam T The type parameter.
 * @param params The parameters to initialize the CalphadBase object.
 */
template <typename T>
CalphadBase<T>::CalphadBase(const Parameters &params) : CalphadBase(params, false) {
  this->CU_ = std::make_shared<CalphadUtils<T>>();
  this->get_parameters();
}

/**
 * @brief Construct a new CalphadBase object.
 *
 * @tparam T The type parameter.
 * @param params The parameters to initialize the CalphadBase object.
 * @param is_KKS A boolean flag indicating whether to initialize KKS.
 */
template <typename T>
CalphadBase<T>::CalphadBase(const Parameters &params, bool is_KKS)
    : params_(params), is_KKS_(is_KKS) {
  this->KKS_ = std::make_shared<KKS<T>>();
  this->get_parameters();
}

/**
 * @brief Get Calphad parameters
 * @remark If KKS is enable, associated parameters are returned.
 *
 * @tparam T
 */
template <typename T>
void CalphadBase<T>::get_parameters() {
  this->element_removed_from_ic_ = this->params_.template get_param_value_or_default<std::string>(
      "element_removed_from_ic", CalphadDefaultConstant::element_removed_from_ic);
  if (this->is_KKS_) {
    this->KKS_->get_parameters(*this);
  }
}

/**
 * @brief High-level method for managing equilibrium calculations
 *
 * @tparam T
 * @param dt The time-step of the simulation
 * @param tp_gf The thermodynamic condition for equilibrium calculations
 * @param chemicalsystem The targeted chemical system
 * @param output_system The output variables
 */
template <typename T>
void CalphadBase<T>::global_execute(
    const int dt, const double time_step, const std::vector<T> &tp_gf,
    const std::vector<std::tuple<std::string, std::string>> &chemicalsystem,
    std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>> &output_system,
    const std::vector<std::tuple<std::vector<std::string>, mfem::Vector>> &previous_output_system,
    std::optional<const std::tuple<std::string, T, T>> phase_field_gf,
    std::optional<const std::vector<T>> tp_gf_old,
    std::optional<const std::vector<std::tuple<std::string, std::string, T, T>>> x_gf,
    std::optional<const std::vector<std::tuple<std::string, std::vector<double>>>> coordinates) {
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
  } else {
    // Verify that all required fields are available for KKS execution
    MFEM_VERIFY(phase_field_gf.has_value(),
                "Error: phase_fields_gf is required for KKS execution.");
    MFEM_VERIFY(tp_gf_old.has_value(), "Error: tp_gf_old is required for KKS execution.");
    MFEM_VERIFY(x_gf.has_value(), "Error: x_gf is required for KKS execution.");
    MFEM_VERIFY(coordinates.has_value(), "Error: coordinates is required for KKS execution.");
    // Execute KKS linearization

    this->KKS_->execute_linearization(*this, dt, time_step, tp_gf, *tp_gf_old, *phase_field_gf,
                                      chemicalsystem, *x_gf, *coordinates);
  }
  // Use specific CALPHAD C++ containers to update output_system
  this->update_outputs(dt, nb_nodes, output_system, previous_output_system);
}

/**
 * @brief Clear containers used to store the results of equilibrium calculations
 *
 * @tparam T
 */
template <typename T>
void CalphadBase<T>::clear_containers() {
  // Clear before filling with new results
  this->diffusion_chemical_potentials_.clear();
  this->chemical_potentials_.clear();
  this->mole_fraction_of_phase_.clear();
  this->elem_mole_fraction_by_phase_.clear();
  this->site_fraction_.clear();
  this->energies_of_phases_.clear();
  this->driving_forces_.clear();
  this->heat_capacity_.clear();
  this->nucleus_.clear();
  this->mobilities_.clear();
  this->error_equilibrium_.clear();
  if (this->is_KKS_) {
    this->KKS_->clear_containers();
  }
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
    const int dt, const size_t nb_nodes,
    std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>> &output_system,
    const std::vector<std::tuple<std::vector<std::string>, mfem::Vector>> &previous_output_system) {
  Catch_Time_Section("CalphadBase<T>::update_outputs");

  ////////////////////
  // Update outputs //
  ////////////////////
  T output(nb_nodes);
  auto get_or_default = [&](const auto &map, const auto &key, auto default_value) {
    if (map.contains(key)) {
      return map.at(key);
    } else {
      return default_value;
    }
  };
  int id_output = -1;
  double default_value = 0.;

  for (auto &[output_infos, output_value] : output_system) {
    id_output++;
    const std::string &output_type = output_infos.back();

    // Fill output with the relevant values
    switch (calphad_outputs::from(output_type)) {
      case calphad_outputs::mu: {
        const std::string &output_elem = output_infos[1];
        for (std::size_t i = 0; i < nb_nodes; ++i) {
          default_value = 0.;
          if (this->error_equilibrium_[i] == CalphadDefaultConstant::error_max) {
            default_value = std::get<1>(previous_output_system[id_output])(i);
          }
          output[i] = get_or_default(this->chemical_potentials_, std::make_tuple(i, output_elem),
                                     default_value);
        }
        break;
      }
      case calphad_outputs::dmu: {
        const std::string &output_elem = output_infos[1];
        for (std::size_t i = 0; i < nb_nodes; ++i) {
          default_value = 0.;
          if (this->error_equilibrium_[i] == CalphadDefaultConstant::error_max) {
            default_value = std::get<1>(previous_output_system[id_output])(i);
          }
          output[i] = get_or_default(this->diffusion_chemical_potentials_,
                                     std::make_tuple(i, output_elem), default_value);
        }
        break;
      }
      case calphad_outputs::xph: {
        const std::string &output_phase = output_infos[1];
        for (std::size_t i = 0; i < nb_nodes; ++i) {
          default_value = 0.;
          if (this->error_equilibrium_[i] == CalphadDefaultConstant::error_max) {
            default_value = std::get<1>(previous_output_system[id_output])(i);
          }
          output[i] = get_or_default(this->mole_fraction_of_phase_,
                                     std::make_tuple(i, output_phase), default_value);
        }
        break;
      }
      case calphad_outputs::xp: {
        const std::string &output_elem = output_infos[1];
        const std::string &output_phase = output_infos[2];
        for (std::size_t i = 0; i < nb_nodes; ++i) {
          default_value = 0.;
          if (this->error_equilibrium_[i] == CalphadDefaultConstant::error_max) {
            default_value = std::get<1>(previous_output_system[id_output])(i);
          }
          output[i] = get_or_default(this->elem_mole_fraction_by_phase_,
                                     std::make_tuple(i, output_phase, output_elem), default_value);
        }
        break;
      }
      case calphad_outputs::y: {
        const std::string &output_cons = output_infos[1];
        const int &output_sub = std::stoi(output_infos[2]);
        const std::string &output_phase = output_infos[3];
        for (std::size_t i = 0; i < nb_nodes; ++i) {
          default_value = 0.;
          if (this->error_equilibrium_[i] == CalphadDefaultConstant::error_max) {
            default_value = std::get<1>(previous_output_system[id_output])(i);
          }
          output[i] = get_or_default(this->site_fraction_,
                                     std::make_tuple(i, output_phase, output_cons, output_sub),
                                     default_value);
        }
        break;
      }
      case calphad_outputs::g: {
        const std::string &output_phase = output_infos[1];
        for (std::size_t i = 0; i < nb_nodes; ++i) {
          default_value = 0.;
          if (this->error_equilibrium_[i] == CalphadDefaultConstant::error_max) {
            default_value = std::get<1>(previous_output_system[id_output])(i);
          }
          output[i] = get_or_default(this->energies_of_phases_,
                                     std::make_tuple(i, output_phase, "G"), default_value);
        }
        break;
      }
      case calphad_outputs::gm: {
        const std::string &output_phase = output_infos[1];
        for (std::size_t i = 0; i < nb_nodes; ++i) {
          default_value = 0.;
          if (this->error_equilibrium_[i] == CalphadDefaultConstant::error_max) {
            default_value = std::get<1>(previous_output_system[id_output])(i);
          }
          output[i] = get_or_default(this->energies_of_phases_,
                                     std::make_tuple(i, output_phase, "GM"), default_value);
        }
        break;
      }
      case calphad_outputs::h: {
        const std::string &output_phase = output_infos[1];
        for (std::size_t i = 0; i < nb_nodes; ++i) {
          default_value = 0.;
          if (this->error_equilibrium_[i] == CalphadDefaultConstant::error_max) {
            default_value = std::get<1>(previous_output_system[id_output])(i);
          }
          output[i] = get_or_default(this->energies_of_phases_,
                                     std::make_tuple(i, output_phase, "H"), default_value);
        }
        break;
      }
      case calphad_outputs::hm: {
        const std::string &output_phase = output_infos[1];
        for (std::size_t i = 0; i < nb_nodes; ++i) {
          default_value = 0.;
          if (this->error_equilibrium_[i] == CalphadDefaultConstant::error_max) {
            default_value = std::get<1>(previous_output_system[id_output])(i);
          }
          output[i] = get_or_default(this->energies_of_phases_,
                                     std::make_tuple(i, output_phase, "HM"), default_value);
        }
        break;
      }
      case calphad_outputs::dgm: {
        const std::string &output_phase = output_infos[1];
        for (std::size_t i = 0; i < nb_nodes; ++i) {
          default_value = 0.;
          if (this->error_equilibrium_[i] == CalphadDefaultConstant::error_max) {
            default_value = std::get<1>(previous_output_system[id_output])(i);
          }
          output[i] = get_or_default(this->driving_forces_, std::make_tuple(i, output_phase),
                                     default_value);
        }
        break;
      }
      case calphad_outputs::cp: {
        for (std::size_t i = 0; i < nb_nodes; ++i) {
          if (this->error_equilibrium_[i] == CalphadDefaultConstant::error_max) {
            output[i] = std::get<1>(previous_output_system[id_output])(i);
          } else {
            output[i] = this->heat_capacity_[i];
          }
        }
        break;
      }
      case calphad_outputs::mob: {
        const std::string &output_phase = output_infos[1];
        const std::string &output_elem = output_infos[2];
        for (std::size_t i = 0; i < nb_nodes; ++i) {
          default_value = -std::numeric_limits<double>::max();
          if (this->error_equilibrium_[i] == CalphadDefaultConstant::error_max) {
            default_value = std::get<1>(previous_output_system[id_output])(i);
          }
          output[i] = get_or_default(this->mobilities_,
                                     std::make_tuple(i, output_phase, output_elem), default_value);
        }
        break;
      }
      case calphad_outputs::nucleus: {
        const std::string &output_phase = output_infos[1];
        for (std::size_t i = 0; i < nb_nodes; ++i) {
          default_value = 0.;
          if (this->error_equilibrium_[i] == CalphadDefaultConstant::error_max) {
            default_value = std::get<1>(previous_output_system[id_output])(i);
          }
          output[i] =
              get_or_default(this->nucleus_, std::make_tuple(i, output_phase), default_value);
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

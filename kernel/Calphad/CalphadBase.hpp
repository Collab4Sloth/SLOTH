
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
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Profiling/Profiling.hpp"

#pragma once

template <typename T>
class CalphadBase {
 protected:
  std::string description_{"UNKNOWN CALPHAD"};

  const Parameters &params_;

  // Containers used to store the results of the equilibrium calculations
  // Chemical potentials for each element. Nodal values
  // key: [node, elem]
  std::map<std::tuple<int, std::string>, double> chemical_potentials_;
  // Element molar fraction for each phase .Nodal values
  // key: [node, phase, elem]
  std::map<std::tuple<int, std::string, std::string>, double> elem_mole_fraction_by_phase_;
  // Site fraction of a given constituant, in a given sublattice for each phase. Nodal values
  // key: [node, phase, cons, sub]
  std::map<std::tuple<int, std::string, std::string, int>, double> site_fraction_;
  // Energy for each phase. Nodal values
  // key: [node, phase, symbol]
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

 public:
  CalphadBase();
  explicit CalphadBase(const Parameters &params);

  virtual void initialize() = 0;

  virtual void execute(const int dt, const std::vector<T> &tp_gf,
                       const std::vector<std::tuple<std::string, std::string>> &chemicalsystem,
                       std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>
                           &output_system) = 0;

  virtual void finalize() = 0;

  ////////////////////////////////

  ////////////////////////////////

  virtual void get_parameters() = 0;

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
CalphadBase<T>::CalphadBase() {}

/**
 * @brief Construct a new Calphad Base< T>:: Calphad Base object
 *
 * @tparam T
 * @param params
 */
template <typename T>
CalphadBase<T>::CalphadBase(const Parameters &params) : params_(params) {}

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
  Catch_Time_Section("Calphad::update_outputs");

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

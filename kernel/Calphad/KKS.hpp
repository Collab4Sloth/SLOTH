
/**
 * @file KKS.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief All needed for Calphad calculation including KKS model
 * @version 0.1
 * @date 2025-05-23
 *
 * @copyright Copyright (c) 2025
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

template <typename T>
class KKS {
 private:
  const double xmin_tol_ = 1.e-12;
  const double kks_abs_tol_solver_ = 1.e-16;
  const double kks_rel_tol_solver_ = 1.e-16;
  const int kks_max_iter_solver_ = 100;
  const int kks_print_level_solver_ = 0;
  std::string KKS_nucleation_strategy_;
  double given_melting_temperature_{std::numeric_limits<double>::max()};
  std::set<int> check_nucleation(CalphadBase<T> &CALPHAD, const std::set<int> &indices_ph_1,
                                 const T &tp_gf_ph_1);

 protected:
  std::string element_removed_from_ic_;
  std::shared_ptr<CalphadUtils<T>> CU_;
  std::string KKS_secondary_phase_;
  double KKS_mobility_for_seed_;
  double KKS_seed_;
  double KKS_seed_radius_;
  double KKS_temperature_increment_;
  double KKS_composition_increment_;
  double KKS_threshold_;
  double KKS_temperature_threshold_;
  bool KKS_temperature_scheme_;
  bool KKS_nucleation_started_{false};
  bool KKS_freeze_nucleation_{true};

  PotentialFunctions<0, ThermodynamicsPotentialDiscretization::Implicit,
                     ThermodynamicsPotentials::H>
      interpolation_func_;

  std::unique_ptr<mfem::SparseMatrix> get_A4linearKKS(
      const std::vector<std::tuple<std::string, std::string>> &chemicalsystem,
      const std::string &phase, const int node);

  mfem::Vector get_h4linearKKS(
      const std::vector<std::tuple<std::string, std::string>> &chemicalsystem,
      const std::string &phase, const int node);

  mfem::Vector get_m4linearKKS(
      const std::vector<std::tuple<std::string, std::string>> &chemicalsystem,
      const std::string &phase, const int node);

  std::map<std::tuple<int, std::string, std::string>, double> chemical_potentials_by_phase_;
  std::map<std::tuple<int, std::string, std::string>, double> chemical_potentials_left_T_;
  std::map<std::tuple<int, int, std::string, std::string>, double> chemical_potentials_left_x_;
  std::map<std::tuple<int, std::string, std::string>, double> chemical_potentials_right_T_;
  std::map<std::tuple<int, int, std::string, std::string>, double> chemical_potentials_right_x_;

 public:
  KKS();
  void get_parameters(const CalphadBase<T> &CALPHAD);
  void execute_linearization(
      CalphadBase<T> &CALPHAD, const int dt, const double time_step, const std::vector<T> &tp_gf,
      const std::vector<T> &tp_gf_old, const std::tuple<std::string, T, T> &phasefields_gf,
      const std::vector<std::tuple<std::string, std::string>> &chemicalsystem,
      const std::vector<std::tuple<std::string, std::string, T, T>> &x_gf,
      const std::vector<std::tuple<std::string, T>> &coord_gf);

  void clear_containers();

  virtual ~KKS();
};
////////////////////////////////
////////////////////////////////

/**
 * @brief Construct a new Calphad Base< T>:: Calphad Base object
 *
 * @tparam T
 */
template <typename T>
KKS<T>::KKS() {
  this->CU_ = std::make_shared<CalphadUtils<T>>();
}

/**
 * @brief Get all parameters required by KKS problem
 *
 * @tparam T
 */
template <typename T>
void KKS<T>::get_parameters(const CalphadBase<T> &CALPHAD) {
  this->element_removed_from_ic_ = CALPHAD.element_removed_from_ic_;

  this->KKS_secondary_phase_ =
      CALPHAD.params_.template get_param_value<std::string>("KKS_secondary_phase");
  this->KKS_temperature_increment_ =
      CALPHAD.params_.template get_param_value<double>("KKS_temperature_increment");
  this->KKS_composition_increment_ =
      CALPHAD.params_.template get_param_value<double>("KKS_composition_increment");
  this->KKS_temperature_scheme_ = CALPHAD.params_.template get_param_value_or_default<bool>(
      "KKS_temperature_explicit_scheme", false);

  this->KKS_seed_ = CALPHAD.params_.template get_param_value_or_default<double>("KKS_seed", 1.);
  this->KKS_seed_radius_ =
      CALPHAD.params_.template get_param_value_or_default<double>("KKS_seed_radius", 1.e-4);
  this->KKS_threshold_ =
      CALPHAD.params_.template get_param_value_or_default<double>("KKS_threshold", 1.e-3);

  this->KKS_temperature_threshold_ = CALPHAD.params_.template get_param_value_or_default<double>(
      "KKS_temperature_threshold", 2800);

  // Nucleation strategy
  this->KKS_nucleation_strategy_ = CALPHAD.params_.template get_param_value_or_default<std::string>(
      "KKS_nucleation_strategy", "LiquidFraction");
  if (KKS_nucleation_strategy::from(this->KKS_nucleation_strategy_) ==
      KKS_nucleation_strategy::given_melting_temperature) {
    this->given_melting_temperature_ =
        CALPHAD.params_.template get_param_value<double>("KKS_given_melting_temperature");
  }

  this->KKS_freeze_nucleation_ =
      CALPHAD.params_.template get_param_value_or_default<bool>("KKS_freeze_nucleation", true);

  this->KKS_nucleation_started_ =
      CALPHAD.params_.template get_param_value_or_default<bool>("KKS_nucleation_started", false);

  this->KKS_mobility_for_seed_ =
      CALPHAD.params_.template get_param_value_or_default<double>("KKS_mobility", 1.);
}

/**
 * @brief Check the nucleation state
 *
 * @tparam T
 * @param indices_ph_1
 * @param tp_gf_ph_1
 * @return std::set<int>
 */
template <typename T>
std::set<int> KKS<T>::check_nucleation(CalphadBase<T> &CALPHAD, const std::set<int> &indices_ph_1,
                                       const T &tp_gf_ph_1) {
  std::set<int> indices_nucleation;

  switch (KKS_nucleation_strategy::from(this->KKS_nucleation_strategy_)) {
    //
    case KKS_nucleation_strategy::liquid_fraction: {
      for (const auto &node : indices_ph_1) {
        // Check only if equilibrium is found
        if (CALPHAD.error_equilibrium_[node] == CalphadDefaultConstant::error_max) continue;
        // Check if secondary phase is found
        if (CALPHAD.elem_mole_fraction_by_phase_.contains(std::make_tuple(
                node, this->KKS_secondary_phase_, this->element_removed_from_ic_))) {
          indices_nucleation.insert(node);
        }
      }
      break;
    }
    //
    case KKS_nucleation_strategy::given_melting_temperature: {
      break;
      for (const auto &node : indices_ph_1) {
        // Check if temperature at node is greater than a given limit
        if (tp_gf_ph_1(node) > this->given_melting_temperature_) {
          indices_nucleation.insert(node);
        }
      }
    }
    default: {
      throw std::runtime_error(
          "KKS<T>::check_nucleation_at_node: error in the choice of method used to compute "
          "the nucleation. Available choices are : liquid_fraction (Default), "
          "given_melting_temperature");
    }
  }
  return indices_nucleation;
}

/**
 * @brief Compute Calphad contributions needed by KKS resolution
 *
 * @tparam T The type of the template parameter (mfem::Vector or std::vector)
 * @param CALPHAD The CALPHAD object for thermodynamic calculation
 * @param dt The time step for the simulation.
 * @param tp_gf The thermodynamic conditions at the current time step.
 * @param tp_gf_old The thermodynamic conditions at the previous time step.
 * @param phasefields_gf The phase-field at the current and previous time-step
 * @param chemicalsystem The targeted chemical system
 * @param x_gf The mole fraction of compoenent at the current and previous time-step
 */
template <typename T>
void KKS<T>::execute_linearization(
    CalphadBase<T> &CALPHAD, const int dt, const double time_step, const std::vector<T> &tp_gf,
    const std::vector<T> &tp_gf_old, const std::tuple<std::string, T, T> &phasefields_gf,
    const std::vector<std::tuple<std::string, std::string>> &chemicalsystem,
    const std::vector<std::tuple<std::string, std::string, T, T>> &x_gf,
    const std::vector<std::tuple<std::string, T>> &coordinates) {
  const int nb_elem = chemicalsystem.size();
  // Creation initial list of nodes
  const size_t nb_nodes = this->CU_->get_size(tp_gf[0]);

  ////////////////////////////////////////////////////////
  // Interpolation: H(phi, t+dt)=Hphi , H(phi, t)=Hphi_old
  ////////////////////////////////////////////////////////
  T Hphi(nb_nodes);
  T Hphi_old(nb_nodes);

  ////////////////////////////////////////////////////////
  // List of nodes by phases and within interface
  ////////////////////////////////////////////////////////
  const auto &[phase, phi_gf, phi_gf_old] = phasefields_gf;

  std::set<int> indices_ph_1;
  std::set<int> indices_ph_2;
  std::set<int> indices_inter;
  for (int i = 0; i < nb_nodes; ++i) {
    const double phi = phi_gf[i];
    if (phi > 1 - this->KKS_threshold_) {
      indices_ph_1.insert(i);
    } else if (phi < this->KKS_threshold_) {
      SlothInfo::debug("Execute_linearization: secondary_phase node ", i);
      indices_ph_2.insert(i);

    } else {
      SlothInfo::debug("Execute_linearization: interfacial node  ", i);
      indices_inter.insert(i);
    }

    // Interpolation function at node
    FType H = this->interpolation_func_.getPotentialFunction(phi_gf_old(i));
    Hphi(i) = H(phi);
    Hphi_old(i) = H(phi_gf_old(i));
  }

  // Lambda for equilibrium calculations performed in the pure phase or in the interface with a
  // temperature deviation
  auto calculate_interface =
      [&](const std::set<int> &indices, const std::vector<T> &delta_tp_gf, const double increment,
          const int id_incr,
          std::map<std::tuple<int, std::string, std::string>, double> &chemical_potential_interface,
          std::vector<std::tuple<std::string, std::string, double>> given_phase,
          const std::string &primary_phase) {
        std::vector<T> delta_tp = delta_tp_gf;
        if (id_incr > -1) {
          delta_tp[id_incr] += increment;
        }
        CALPHAD.execute(dt, indices, delta_tp, chemicalsystem, given_phase);

        for (const auto &in : indices) {
          if (CALPHAD.error_equilibrium_[in] == CalphadDefaultConstant::error_max) continue;
          for (const auto &elem : chemicalsystem) {
            const auto &[elem1, unit] = elem;
            const double mu = CALPHAD.chemical_potentials_.at(std::make_tuple(in, elem1));
            chemical_potential_interface.emplace(std::make_tuple(in, elem1, primary_phase), mu);
          }
        }
      };
  //
  // Lambda for equilibrium calculations performed in the pure phase with a deviation in
  // composition
  auto calculate_interface_x =
      [&](const std::vector<T> &delta_tp_gf, const double increment, const int id_incr,
          int index_el,
          std::map<std::tuple<int, int, std::string, std::string>, double>
              &chemical_potential_interface,
          std::vector<std::tuple<std::string, std::string, double>> given_phase,
          const std::string &primary_phase) {
        std::vector<T> delta_tp = delta_tp_gf;
        if (id_incr > -1) {
          delta_tp[id_incr] += increment;
        }

        CALPHAD.execute(dt, indices_inter, delta_tp, chemicalsystem, given_phase);
        for (const auto &in : indices_inter) {
          if (CALPHAD.error_equilibrium_[in] == CalphadDefaultConstant::error_max) continue;
          for (const auto &elem : chemicalsystem) {
            const auto &[elem1, unit] = elem;
            const double mu = CALPHAD.chemical_potentials_.at(std::make_tuple(in, elem1));
            chemical_potential_interface.emplace(
                std::make_tuple(index_el, in, elem1, primary_phase), mu);
          }
        }
      };

  ////////////////////////////////////////////////////////
  /// Point of expansion (Tbar,Xbar)
  ////////////////////////////////////////////////////////

  // Primary Phase : index _1
  std::vector<T> bar_tp_gf_ph_1(tp_gf.size());
  std::vector<T> pure_bar_tp_gf_ph_1(tp_gf.size());
  // KKS_secondary_phase : index _2
  std::vector<T> bar_tp_gf_ph_2(tp_gf.size());
  std::vector<T> pure_bar_tp_gf_ph_2(tp_gf.size());
  // Tbar
  if (!this->KKS_temperature_scheme_) {
    bar_tp_gf_ph_1[0] = tp_gf[0];
    bar_tp_gf_ph_2[0] = tp_gf[0];
    pure_bar_tp_gf_ph_1[0] = tp_gf[0];
    pure_bar_tp_gf_ph_2[0] = tp_gf[0];
  } else {
    bar_tp_gf_ph_1[0] = tp_gf_old[0];
    bar_tp_gf_ph_2[0] = tp_gf_old[0];
    pure_bar_tp_gf_ph_1[0] = tp_gf_old[0];
    pure_bar_tp_gf_ph_2[0] = tp_gf_old[0];
  }
  // Pressure
  bar_tp_gf_ph_1[1] = tp_gf[1];
  bar_tp_gf_ph_2[1] = tp_gf[1];
  pure_bar_tp_gf_ph_1[1] = tp_gf[1];
  pure_bar_tp_gf_ph_2[1] = tp_gf[1];

  // Composition : Xbar
  MFEM_VERIFY(x_gf.size() > 0, "Error while getting molar fraction by phase");

  for (const auto &[elem, phase_elem, x_elem, x_elem_old] : x_gf) {
    auto it = std::ranges::find_if(chemicalsystem,
                                   [&elem](const auto &t) { return std::get<0>(t) == elem; });

    std::optional<int> id;
    if (it != chemicalsystem.end()) {
      id = std::distance(chemicalsystem.begin(), it);
    }
    MFEM_VERIFY(id.has_value(), "Error while getting the element");
    const int id_value = id.value();
    const bool is_primary_phase = (phase_elem == phase);
    const bool is_secondary_phase = (phase_elem == this->KKS_secondary_phase_);

    if (is_primary_phase) {
      if (elem != this->element_removed_from_ic_) {
        bar_tp_gf_ph_1[id_value + 2] = x_elem;
      } else {
        bar_tp_gf_ph_1[id_value + 2] = tp_gf[id_value + 2];
      }
      pure_bar_tp_gf_ph_1[id_value + 2] = tp_gf[id_value + 2];

    } else if (is_secondary_phase) {
      if (elem != this->element_removed_from_ic_) {
        bar_tp_gf_ph_2[id_value + 2] = x_elem;

        // When scecondary phase is previously absent, initialize with primary phase
        if (!indices_inter.empty()) {
          for (const int id_new : indices_inter) {
            if (bar_tp_gf_ph_2[id_value + 2](id_new) < this->xmin_tol_) {
              bar_tp_gf_ph_2[id_value + 2](id_new) = tp_gf[id_value + 2](id_new);
            }
          }
        }
      } else {
        bar_tp_gf_ph_2[id_value + 2] = tp_gf[id_value + 2];
      }
      pure_bar_tp_gf_ph_2[id_value + 2] = tp_gf[id_value + 2];

    } else {
      const std::string &error_msg =
          "Error while setting molar fraction at point approximation. Unknown phase " + phase_elem;
      mfem::mfem_error(error_msg.c_str());
    }
  }

  ////////////////////////////////////////////////////////
  /// Calculation of all thermodynamic contributions
  ////////////////////////////////////////////////////////

  /////////////////////////
  /// Primary phase
  /////////////////////////
  std::vector<std::tuple<std::string, std::string, double>> st_phase_12 = {
      {phase, "entered", 0.}, {this->KKS_secondary_phase_, "entered", 0.}};

  if ((this->KKS_nucleation_started_ && this->KKS_freeze_nucleation_) ||
      (pure_bar_tp_gf_ph_1[0][0] < this->KKS_temperature_threshold_)) {
    st_phase_12 = {{phase, "entered", 0.}};
  }

  // In practice, two entered phase seems to be more stable than considering dormant phase
  // {phase, "entered"}, {this->KKS_secondary_phase_, "dormant"}};
  std::vector<std::tuple<std::string, std::string, double>> st_phase_1 = {{phase, "entered", 0.}};
  std::vector<std::tuple<std::string, std::string, double>> st_phase_2 = {
      {this->KKS_secondary_phase_, "entered", 0.}};

  calculate_interface(indices_ph_1, pure_bar_tp_gf_ph_1, 0., -1,
                      this->chemical_potentials_by_phase_, st_phase_12, phase);

  // Cancel before checking nucleation state
  for (int i = 0; i < nb_nodes; ++i) {
    // Must be zero except for nodes detected as  nucleus  in solid
    CALPHAD.nucleus_[std::make_tuple(i, this->KKS_secondary_phase_)] = 0.;
    // Must be zero except for nodes in interface.
    CALPHAD.driving_forces_[std::make_tuple(i, this->KKS_secondary_phase_)] = 0.;
  }

  if (!this->KKS_nucleation_started_) {
    bool local_nucleation = false;
    std::set<int> indices_nucleation =
        this->check_nucleation(CALPHAD, indices_ph_1, pure_bar_tp_gf_ph_1[0]);

    // Nucleation detected ?
    if (indices_nucleation.size() > 0) {
      local_nucleation = true;
      // If nucleation is not frozen once detected, local_nucleation must be false
      local_nucleation &= this->KKS_freeze_nucleation_;
    }

    // Create circular nucleus around the node where secondary phase is detected
    // No need to check error_equilibrium because of already done when indices_nucleation is created
    for (const auto &inuc : indices_nucleation) {
      SlothInfo::debug("Nucleation detected at node ", inuc, " T ", pure_bar_tp_gf_ph_1[0](inuc));
      // const double seed =
      //     CALPHAD.mole_fraction_of_phase_[std::make_tuple(inuc, this->KKS_secondary_phase_)];
      const double seed = this->KKS_seed_;
      const double nucleus_value = -seed / (time_step * this->KKS_mobility_for_seed_);
      for (const auto &j : indices_ph_1) {
        double rr = 0;

        for (const auto &[_, coord] : coordinates) {
          const double coord_diff = coord[inuc] - coord[j];
          rr += coord_diff * coord_diff;
        }
        if (std::sqrt(rr) < this->KKS_seed_radius_) {
          SlothInfo::debug("Nucleation extended around node ", inuc, " with node ", j);
          // double seed =
          //     CALPHAD.mole_fraction_of_phase_[std::make_tuple(node, this->KKS_secondary_phase_)];
          // To enhance the chance to initiate the phase change
          CALPHAD.nucleus_[std::make_tuple(j, this->KKS_secondary_phase_)] = nucleus_value;
        }
      }
    }

    // Synchronization of nucleation_started across all MPI ranks
    int local_flag = (local_nucleation ? 1 : 0);
    int global_flag = 0;
    MPI_Allreduce(&local_flag, &global_flag, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    this->KKS_nucleation_started_ = (global_flag == 1);
  }

  /////////////////////////
  /// Secondary phase (if exists)
  /////////////////////////
  if (!indices_ph_2.empty()) {
    calculate_interface(indices_ph_2, pure_bar_tp_gf_ph_2, 0., -1,
                        this->chemical_potentials_by_phase_, st_phase_2,
                        this->KKS_secondary_phase_);
  }

  /////////////////////////
  /// Interface (if exists)
  /////////////////////////
  // 6 + 4 * (nelem-1) times more equilibrium calculation
  //
  if (!indices_inter.empty()) {
    // T : equilibrium calculations in both phases (Tbar, Xbar)
    calculate_interface(indices_inter, bar_tp_gf_ph_1, 0, -1, this->chemical_potentials_by_phase_,
                        st_phase_1, phase);
    calculate_interface(indices_inter, bar_tp_gf_ph_2, 0., -1, this->chemical_potentials_by_phase_,
                        st_phase_2, this->KKS_secondary_phase_);

    // T-dT : equilibrium calculations in both phases (Tbar - deltaT, Xbar)
    calculate_interface(indices_inter, bar_tp_gf_ph_1, -this->KKS_temperature_increment_, 0,
                        this->chemical_potentials_left_T_, st_phase_1, phase);

    calculate_interface(indices_inter, bar_tp_gf_ph_2, -this->KKS_temperature_increment_, 0,
                        this->chemical_potentials_left_T_, st_phase_2, this->KKS_secondary_phase_);

    // T+dT : equilibrium calculations in both phases (Tbar + deltaT, Xbar)
    calculate_interface(indices_inter, bar_tp_gf_ph_1, this->KKS_temperature_increment_, 0,
                        this->chemical_potentials_right_T_, st_phase_1, phase);

    calculate_interface(indices_inter, bar_tp_gf_ph_2, this->KKS_temperature_increment_, 0,
                        this->chemical_potentials_right_T_, st_phase_2, this->KKS_secondary_phase_);

    // Loop over chemical system (except the reference element)
    int ielem = 0;
    for (const auto &[elem, unit] : chemicalsystem) {
      if (elem != this->element_removed_from_ic_) {
        // x + dx
        // Equilibrium calculations in both phases at (Tbar, Xbar + delta X)
        calculate_interface_x(bar_tp_gf_ph_1, this->KKS_composition_increment_, ielem + 2, ielem,
                              this->chemical_potentials_right_x_, st_phase_1, phase);

        calculate_interface_x(bar_tp_gf_ph_2, this->KKS_composition_increment_, ielem + 2, ielem,
                              this->chemical_potentials_right_x_, st_phase_2,
                              this->KKS_secondary_phase_);

        // x - dx
        // Equilibrium calculations in both phases at (Tbar, Xbar - delta X)
        calculate_interface_x(bar_tp_gf_ph_1, -this->KKS_composition_increment_, ielem + 2, ielem,
                              this->chemical_potentials_left_x_, st_phase_1, phase);

        calculate_interface_x(bar_tp_gf_ph_2, -this->KKS_composition_increment_, ielem + 2, ielem,
                              this->chemical_potentials_left_x_, st_phase_2,
                              this->KKS_secondary_phase_);
      }
      ielem++;
    }

    ////////////////////////////////////////////////////////
    /// Solve linear systems
    ////////////////////////////////////////////////////////
    // Allocate
    mfem::Array<int> offsets(3);
    offsets[0] = 0;
    offsets[1] = nb_elem - 1;
    offsets[2] = 2 * nb_elem - 2;

    mfem::BlockVector bb(offsets);
    mfem::BlockVector deltaX_phase(offsets);
    mfem::Vector HphiNode(nb_elem - 1);
    mfem::Vector one_minus_HphiNode(nb_elem - 1);

    std::unique_ptr<mfem::BlockMatrix> AA = std::make_unique<mfem::BlockMatrix>(offsets, offsets);

    auto solver = std::make_shared<mfem::GMRESSolver>();
    solver->SetOperator(*AA);
    solver->SetAbsTol(this->kks_abs_tol_solver_);
    solver->SetRelTol(this->kks_rel_tol_solver_);
    solver->SetMaxIter(this->kks_max_iter_solver_);
    solver->SetPrintLevel(this->kks_print_level_solver_);

    for (const auto &node : indices_inter) {
      // Solve KKS only if an equilibrium is found. Otherwise, solution at previous time-step will
      // be taken into account
      if (CALPHAD.error_equilibrium_[node] == CalphadDefaultConstant::error_max) continue;

      //
      // Build system
      //
      // Matrix for primary and secondary phase

      std::unique_ptr<mfem::SparseMatrix> Al =
          this->get_A4linearKKS(chemicalsystem, this->KKS_secondary_phase_, node);
      std::unique_ptr<mfem::SparseMatrix> As = this->get_A4linearKKS(chemicalsystem, phase, node);
      // Vector for primary and secondary phase
      mfem::Vector hs = this->get_h4linearKKS(chemicalsystem, phase, node);
      mfem::Vector hl = this->get_h4linearKKS(chemicalsystem, this->KKS_secondary_phase_, node);
      mfem::Vector ms = this->get_m4linearKKS(chemicalsystem, phase, node);
      mfem::Vector ml = this->get_m4linearKKS(chemicalsystem, this->KKS_secondary_phase_, node);

      // ml - ms + (T - Tbar) * (hl - hs)
      // ml - ms
      mfem::Vector ml_minus_ms(ml);
      ml_minus_ms -= ms;

      // T - Tbar
      const double deltaT = tp_gf[0](node) - bar_tp_gf_ph_1[0](node);

      // hl * (T - Tbar)
      mfem::Vector hl_deltaT(hl);
      hl_deltaT *= deltaT;

      // hs * (T - Tbar)
      mfem::Vector hs_deltaT(hs);
      hs_deltaT *= deltaT;

      // hl * (T - Tbar) - hs * (T - Tbar)
      mfem::Vector hl_minus_hs(hl_deltaT);
      hl_minus_hs -= hs_deltaT;

      // LHS : deltaX
      mfem::Vector deltaX(nb_elem - 1);
      int ielem = 0;
      int jelem = 0;
      for (const auto &[elem, unit] : chemicalsystem) {
        if (elem != this->element_removed_from_ic_) {
          const double current_x = tp_gf[jelem + 2](node);
          const double old_x = tp_gf_old[jelem + 2](node);
          const double dx = current_x - old_x;

          const double current_Hphi = Hphi(node);
          const double old_Hphi = Hphi_old(node);

          const double xdH = (bar_tp_gf_ph_1[jelem + 2](node) - bar_tp_gf_ph_2[jelem + 2](node)) *
                             (current_Hphi - old_Hphi);

          SlothInfo::debug("Calculation of dx at node ", node, " = ", dx, " for elem ", elem,
                           " with x ", current_x, " and x_old ", old_x);
          SlothInfo::debug("Calculation of xdH at node ", node, " = ", xdH, " for elem ", elem,
                           " with xs_old ", bar_tp_gf_ph_1[jelem + 2](node), " xl_old ",
                           bar_tp_gf_ph_2[jelem + 2](node), " Hphi(node) ", current_Hphi,
                           " and Hphi_old(node) ", old_Hphi);

          deltaX(ielem) = dx - xdH;

          ielem++;
        }
        jelem++;
      }

      //
      // Assemble system
      //

      // mfem::BlockVector bb(offsets);
      mfem::Vector &b0 = bb.GetBlock(0);
      b0 = ml_minus_ms;
      b0 += hl_minus_hs;

      mfem::Vector &b1 = bb.GetBlock(1);
      b1 = deltaX;
      mfem::Vector &d0 = deltaX_phase.GetBlock(0);
      mfem::Vector &d1 = deltaX_phase.GetBlock(1);
      d0 = 0.;
      d1 = 0.;

      // mfem::BlockMatrix *AA = new mfem::BlockMatrix(offsets, offsets);
      AA->SetBlock(0, 0, As.get());
      *Al *= -1.;
      AA->SetBlock(0, 1, Al.get());

      // identity
      HphiNode = Hphi(node);
      one_minus_HphiNode = 1. - Hphi(node);

      auto A10 = std::make_unique<mfem::SparseMatrix>(HphiNode);
      auto A11 = std::make_unique<mfem::SparseMatrix>(one_minus_HphiNode);

      AA->SetBlock(1, 0, A10.get());
      AA->SetBlock(1, 1, A11.get());

      AA->Finalize();

      //
      // Solve system
      //
      // mfem::BlockVector deltaX_phase(offsets);

      solver->Mult(bb, deltaX_phase);
      // Check result of linear system
      if (Verbosity::Debug <= verbosityLevel) {
        SlothInfo::print("KKS linear system at node (A X = B):", node);
        AA->Print();
        deltaX_phase.Print();
        bb.Print();
      }

      // Back to Al. Required before recovering thermodynamic contribution
      *Al *= -1.;

      ////////////////////////////////////////////////////////
      /// Recover thermodynamic contribution
      ////////////////////////////////////////////////////////
      const mfem::Vector delta_XS = deltaX_phase.GetBlock(0);
      const mfem::Vector delta_XL = deltaX_phase.GetBlock(1);
      ////////////////////////////////////////
      // Composition by phase
      ////////////////////////////////////////
      int i = 0;
      int ie = 0;
      double sumS = 0.;
      double sumL = 0.;
      for (const auto &[elem, unit] : chemicalsystem) {
        if (elem != this->element_removed_from_ic_) {
          CALPHAD.elem_mole_fraction_by_phase_[std::make_tuple(node, phase, elem)] =
              bar_tp_gf_ph_1[i + 2](node) + delta_XS(ie);
          CALPHAD.elem_mole_fraction_by_phase_[std::make_tuple(
              node, this->KKS_secondary_phase_, elem)] = bar_tp_gf_ph_2[i + 2](node) + delta_XL(ie);
          sumS += CALPHAD.elem_mole_fraction_by_phase_[std::make_tuple(node, phase, elem)];
          sumL += CALPHAD.elem_mole_fraction_by_phase_[std::make_tuple(
              node, this->KKS_secondary_phase_, elem)];
          ie++;
        }
        ++i;
      }
      // Last component in each phase
      CALPHAD.elem_mole_fraction_by_phase_[std::make_tuple(
          node, phase, this->element_removed_from_ic_)] = 1. - sumS;
      CALPHAD.elem_mole_fraction_by_phase_[std::make_tuple(
          node, this->KKS_secondary_phase_, this->element_removed_from_ic_)] = 1. - sumL;

      ////////////////////////////////////////
      // Difference of chemical potentials
      ////////////////////////////////////////
      mfem::Vector dmu_s(ms);
      dmu_s += hs_deltaT;

      mfem::Vector AsDeltaXs(nb_elem - 1);
      As->Mult(delta_XS, AsDeltaXs);
      dmu_s += AsDeltaXs;

      mfem::Vector dmu_l(ml);
      dmu_l += hl_deltaT;
      mfem::Vector AlDeltaXl(nb_elem - 1);
      Al->Mult(delta_XL, AlDeltaXl);
      dmu_l += AlDeltaXl;

      ////////////////////////////////////////
      // Gibbs energies (implicit case)
      ////////////////////////////////////////
      double gs = CALPHAD.energies_of_phases_[std::make_tuple(node, phase, "GM")];
      double gl =
          CALPHAD.energies_of_phases_[std::make_tuple(node, this->KKS_secondary_phase_, "GM")];
      auto DeltaXsAsDeltaXs = AsDeltaXs * delta_XS;
      DeltaXsAsDeltaXs *= 0.5;
      auto DeltaXlAlDeltaXl = AlDeltaXl * delta_XL;
      DeltaXlAlDeltaXl *= 0.5;
      auto DeltaXsMs = ms * delta_XS;
      auto DeltaXlMl = ml * delta_XL;
      gs += DeltaXsMs;
      gs += DeltaXsAsDeltaXs;
      gl += DeltaXlMl;
      gl += DeltaXlAlDeltaXl;

      // Molar Gibbs Energy by phase
      CALPHAD.energies_of_phases_[std::make_tuple(node, phase, "GM")] = gs;
      CALPHAD.energies_of_phases_[std::make_tuple(node, this->KKS_secondary_phase_, "GM")] = gl;

      // Driving force by phase : contribution -gm
      CALPHAD.driving_forces_[std::make_tuple(node, phase)] =
          -CALPHAD.energies_of_phases_[std::make_tuple(node, phase, "GM")];
      CALPHAD.driving_forces_[std::make_tuple(node, this->KKS_secondary_phase_)] =
          -CALPHAD.energies_of_phases_[std::make_tuple(node, this->KKS_secondary_phase_, "GM")];

      i = 0;
      ie = 0;
      for (const auto &[elem, unit] : chemicalsystem) {
        if (elem != this->element_removed_from_ic_) {
          SlothInfo::debug("Check at node ", node, " for elem ", elem);
          SlothInfo::debug("deltaX to obtain ", deltaX(ie), " deltaX predicted ",
                           delta_XS(ie) * Hphi(node) + delta_XL(ie) * (1. - Hphi(node)));
          SlothInfo::debug(
              "x to obtain ", tp_gf[i + 2](node), " x predicted ",
              CALPHAD.elem_mole_fraction_by_phase_[std::make_tuple(node, phase, elem)] *
                      Hphi(node) +
                  CALPHAD.elem_mole_fraction_by_phase_[std::make_tuple(
                      node, this->KKS_secondary_phase_, elem)] *
                      (1. - Hphi(node)));

          SlothInfo::debug("Diffusion chemical potential: primary phase ", dmu_s(ie),
                           " ; secondary phase", dmu_l(ie));

          // Diffusion chemical potential
          CALPHAD.diffusion_chemical_potentials_[std::make_tuple(node, elem)] = dmu_s(ie);

          // Pseudo Driving force : contribution dmu_i * x_i
          //
          CALPHAD.driving_forces_[std::make_tuple(node, phase)] +=
              CALPHAD.diffusion_chemical_potentials_[std::make_tuple(node, elem)] *
              CALPHAD.elem_mole_fraction_by_phase_[std::make_tuple(node, phase, elem)];

          CALPHAD.driving_forces_[std::make_tuple(node, this->KKS_secondary_phase_)] +=
              CALPHAD.diffusion_chemical_potentials_[std::make_tuple(node, elem)] *
              CALPHAD.elem_mole_fraction_by_phase_[std::make_tuple(node, this->KKS_secondary_phase_,
                                                                   elem)];
          ie++;
        }
        i++;
      }

      SlothInfo::debug("Gibbs energy: primary phase ", gs, " ; secondary phase ", gl);
    }
  }
}

/**
 * @brief Second order finite difference approximation
 *
 * @tparam T The type of the template parameter (mfem::Vector or std::vector)
 * @param chemicalsystem The targeted chemical system
 * @param phase The phase for which derivatives are computed
 * @param node The node where calculation are done
 * @return mfem::SparseMatrix*
 */
template <typename T>
std::unique_ptr<mfem::SparseMatrix> KKS<T>::get_A4linearKKS(
    const std::vector<std::tuple<std::string, std::string>> &chemicalsystem,
    const std::string &phase, const int node) {
  const int nb_elem = chemicalsystem.size();
  auto AA = std::make_unique<mfem::SparseMatrix>(nb_elem - 1, nb_elem - 1);

  //  Diagonal
  // d2dxdx =  d(mu_x - mu-n)/dx  =[ (mu_x - mu-n)(x+dx) - (mu_x - mu-n)(x-dx) ] / 2dx
  // elem_id corresponds to index of variable with delta_x, index of derivative
  int elem_id = 0;
  int vid = 0;
  for (const auto &[elem, unit] : chemicalsystem) {
    if (elem != this->element_removed_from_ic_) {
      double d2gd2x =
          this->chemical_potentials_right_x_[std::make_tuple(elem_id, node, elem, phase)] -
          this->chemical_potentials_right_x_[std::make_tuple(
              elem_id, node, this->element_removed_from_ic_, phase)] -
          (this->chemical_potentials_left_x_[std::make_tuple(elem_id, node, elem, phase)] -
           this->chemical_potentials_left_x_[std::make_tuple(
               elem_id, node, this->element_removed_from_ic_, phase)]);
      d2gd2x /= 2. * this->KKS_composition_increment_;
      AA->Set(vid, vid, d2gd2x);
      vid++;
    }
    elem_id++;
  }
  // Off-diagonal (and symmetry)
  elem_id = 0;
  vid = 0;
  for (const auto &[ielem, unit] : chemicalsystem) {
    if (ielem != this->element_removed_from_ic_) {
      int elem_jd = 0;
      int vjd = 0;
      ///////////
      // row i
      for (const auto &[jelem, unit] : chemicalsystem) {
        if (jelem != this->element_removed_from_ic_) {
          ///////////
          // col j
          //
          // d2dxdy = 0.5 d(mu_x - mu-n)/dy + 0.5 d(mu_y - mu-n)/dx
          if (ielem != jelem) {
            double dmuxdy =
                (this->chemical_potentials_right_x_[std::make_tuple(elem_jd, node, ielem, phase)] -
                 this->chemical_potentials_right_x_[std::make_tuple(
                     elem_jd, node, this->element_removed_from_ic_, phase)]) -
                ((this->chemical_potentials_left_x_[std::make_tuple(elem_jd, node, ielem, phase)] -
                  this->chemical_potentials_left_x_[std::make_tuple(
                      elem_jd, node, this->element_removed_from_ic_, phase)]));
            dmuxdy /= 2. * this->KKS_composition_increment_;
            double dmuydx =
                (this->chemical_potentials_right_x_[std::make_tuple(elem_id, node, jelem, phase)] -
                 this->chemical_potentials_right_x_[std::make_tuple(
                     elem_id, node, this->element_removed_from_ic_, phase)]) -
                ((this->chemical_potentials_left_x_[std::make_tuple(elem_id, node, jelem, phase)] -
                  this->chemical_potentials_left_x_[std::make_tuple(
                      elem_id, node, this->element_removed_from_ic_, phase)]));
            dmuydx /= 2. * this->KKS_composition_increment_;

            double d2gdxdy = 0.5 * dmuxdy + 0.5 * dmuydx;

            AA->Set(vid, vjd, d2gdxdy);
          }
          vjd++;
          ///////////
        }
        elem_jd++;
      }
      ///////////
      vid++;
    }
    elem_id++;
  }

  return AA;
}

/**
 * @brief First order finite difference approximation
 *
 * @tparam T The type of the template parameter (mfem::Vector or std::vector)
 * @param chemicalsystem The targeted chemical system
 * @param phase The phase for which derivatives are computed
 * @param node The node where calculation are done
 * @return mfem::Vector
 */
template <typename T>
mfem::Vector KKS<T>::get_m4linearKKS(
    const std::vector<std::tuple<std::string, std::string>> &chemicalsystem,
    const std::string &phase, const int node) {
  const int nb_elem = chemicalsystem.size();
  mfem::Vector mm(nb_elem - 1);
  int vid = 0;

  for (const auto &[ielem, unit] : chemicalsystem) {
    if (ielem != this->element_removed_from_ic_) {
      mm(vid) = this->chemical_potentials_by_phase_[std::make_tuple(node, ielem, phase)] -
                this->chemical_potentials_by_phase_[std::make_tuple(
                    node, this->element_removed_from_ic_, phase)];
      vid++;
    }
  }

  return mm;
}

/**
 * @brief Second order finite difference approximation (cross derivatives)
 *
 * @tparam T The type of the template parameter (mfem::Vector or std::vector)
 * @param chemicalsystem The targeted chemical system
 * @param phase The phase for which derivatives are computed
 * @param node The node where calculation are done
 * @return mfem::Vector
 */
template <typename T>
mfem::Vector KKS<T>::get_h4linearKKS(
    const std::vector<std::tuple<std::string, std::string>> &chemicalsystem,
    const std::string &phase, const int node) {
  const int nb_elem = chemicalsystem.size();
  mfem::Vector hh(nb_elem - 1);
  int vid = 0;
  for (const auto &[ielem, unit] : chemicalsystem) {
    if (ielem != this->element_removed_from_ic_) {
      hh(vid) = this->chemical_potentials_right_T_[std::make_tuple(node, ielem, phase)] -
                this->chemical_potentials_right_T_[std::make_tuple(
                    node, this->element_removed_from_ic_, phase)] -
                (this->chemical_potentials_left_T_[std::make_tuple(node, ielem, phase)] -
                 this->chemical_potentials_left_T_[std::make_tuple(
                     node, this->element_removed_from_ic_, phase)]);
      hh(vid) /= 2.0 * this->KKS_temperature_increment_;
      vid++;
    }
  }

  return hh;
}

/**
 * @brief Clear containers used to store the results of equilibrium calculations
 *
 * @tparam T
 */
template <typename T>
void KKS<T>::clear_containers() {
  // Clear before filling with new results
  this->chemical_potentials_by_phase_.clear();
  this->chemical_potentials_left_T_.clear();
  this->chemical_potentials_left_x_.clear();
  this->chemical_potentials_right_T_.clear();
  this->chemical_potentials_right_x_.clear();
}

/**
 * @brief Destroy the Calphad Base< T>:: Calphad Base object
 *
 * @tparam T
 */
template <typename T>
KKS<T>::~KKS() {}

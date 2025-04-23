
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

// TODO(ci) ajouter des conditions kks pour etre obligatoirement en x, avoir acces au dernier
// element, à la nucleation, à delta_x, delta_t

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
  std::string element_removed_from_ic_;
  const Parameters &params_;

  // KKS
  mfem::SparseMatrix get_A4linearKKS(
      const std::vector<std::tuple<std::string, std::string>> &chemicalsystem,
      const std::string &phase, const int node);
  mfem::Vector get_h4linearKKS(
      const std::vector<std::tuple<std::string, std::string>> &chemicalsystem,
      const std::string &phase, const int node);
  mfem::Vector get_m4linearKKS(
      const std::vector<std::tuple<std::string, std::string>> &chemicalsystem,
      const std::string &phase, const int node);

  // Containers used to store the results of the equilibrium calculations
  // Chemical potentials for each element. Nodal values

  // Chemical potential. Nodal values
  // key: [node, elem]
  std::map<std::tuple<int, std::string>, double> chemical_potentials_;
  std::map<std::tuple<int, std::string, std::string>, double> chemical_potentials_by_phase_;
  std::map<std::tuple<int, std::string, std::string>, double> chemical_potentials_left_T_;
  std::map<std::tuple<int, int, std::string, std::string>, double> chemical_potentials_left_x_;
  std::map<std::tuple<int, std::string, std::string>, double> chemical_potentials_right_T_;
  std::map<std::tuple<int, int, std::string, std::string>, double> chemical_potentials_right_x_;
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
  this->element_removed_from_ic_ = this->params_.template get_param_value_or_default<std::string>(
      "element_removed_from_ic", CalphadDefaultConstant::element_removed_from_ic);
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

  //
  MFEM_VERIFY(this->params_.has_parameter("element_removed_from_ic"),
              "Error while defining Calphad object with KKS.  Parameter element_removed_from_ic "
              "must be defined. Please check your data.");
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
  // TODO(cci) add nucleation managment
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

  // Interface calculations
  auto calculate_interface = [&](const std::set<int> &indices, const std::vector<T> &delta_tp_gf,
                                 const std::string &given_phase,
                                 std::map<std::tuple<int, std::string, std::string>, double>
                                     &chemical_potential_interface) {
    this->execute(dt, indices, delta_tp_gf, chemicalsystem, given_phase);
    for (const auto &in : indices) {
      for (const auto &elem : chemicalsystem) {
        const auto &[elem1, unit] = elem;
        const double mu = this->chemical_potentials_.at(std::make_tuple(in, elem1));
        chemical_potential_interface.emplace(std::make_tuple(in, elem1, given_phase), mu);
      }
    }
  };
  //
  auto calculate_interface_x = [&](const std::vector<T> &delta_tp_gf,
                                   const std::string &given_phase, const int index_el,
                                   std::map<std::tuple<int, int, std::string, std::string>, double>
                                       &chemical_potential_interface) {
    this->execute(dt, indices_inter, delta_tp_gf, chemicalsystem, given_phase);
    for (const auto &in : indices_inter) {
      for (const auto &elem : chemicalsystem) {
        const auto &[elem1, unit] = elem;
        const double mu = this->chemical_potentials_.at(std::make_tuple(in, elem1));
        chemical_potential_interface.emplace(std::make_tuple(index_el, in, elem1, given_phase), mu);
      }
    }
  };
  // T
  calculate_interface(indices_ph_1, tp_gf, phase, this->chemical_potentials_by_phase_);
  calculate_interface(indices_ph_2, tp_gf, this->KKS_secondary_phase_,
                      this->chemical_potentials_by_phase_);

  // T-dT
  std::vector<T> delta_tp_gf = tp_gf;
  delta_tp_gf[0] += this->KKS_temperature_increment_;
  calculate_interface(indices_inter, delta_tp_gf, phase, this->chemical_potentials_left_T_);
  calculate_interface(indices_inter, delta_tp_gf, this->KKS_secondary_phase_,
                      this->chemical_potentials_left_T_);

  // T+dT
  delta_tp_gf = tp_gf;
  delta_tp_gf[0] -= this->KKS_temperature_increment_;
  calculate_interface(indices_inter, delta_tp_gf, phase, this->chemical_potentials_right_T_);
  calculate_interface(indices_inter, delta_tp_gf, this->KKS_secondary_phase_,
                      this->chemical_potentials_right_T_);

  // Loop over chemicalsystem
  for (std::size_t i = 0; i < chemicalsystem.size(); ++i) {
    // x+dx
    delta_tp_gf = tp_gf;
    delta_tp_gf[i + 2] += this->KKS_composition_increment_;
    calculate_interface_x(delta_tp_gf, phase, i, this->chemical_potentials_right_x_);
    calculate_interface_x(delta_tp_gf, this->KKS_secondary_phase_, i,
                          this->chemical_potentials_right_x_);

    // x-dx
    delta_tp_gf = tp_gf;
    delta_tp_gf[i + 2] -= this->KKS_composition_increment_;
    calculate_interface_x(delta_tp_gf, phase, i, this->chemical_potentials_left_x_);
    calculate_interface_x(delta_tp_gf, this->KKS_secondary_phase_, i,
                          this->chemical_potentials_left_x_);
  }
  ////////////////////////////////////////////////////////
  /// Solve linear systems
  ////////////////////////////////////////////////////////
  auto ss = std::make_shared<mfem::HypreGMRES>(MPI_COMM_WORLD);

  int node = 0;

  mfem::SparseMatrix Al = get_A4linearKKS(chemicalsystem, this->KKS_secondary_phase_, node);
  mfem::SparseMatrix As = get_A4linearKKS(chemicalsystem, phase, node);
  mfem::Vector hs = get_h4linearKKS(chemicalsystem, phase, node);
  mfem::Vector hl = get_h4linearKKS(chemicalsystem, this->KKS_secondary_phase_, node);
  mfem::Vector ms = get_m4linearKKS(chemicalsystem, phase, node);
  mfem::Vector ml = get_m4linearKKS(chemicalsystem, this->KKS_secondary_phase_, node);
  mfem::Vector b;
  b = ms;
  b -= ml;

  // mfem::SparseMatrix A(3);

  // A.Add(i, j, mu);

  // A.Finalize();

  // mfem::HypreParVector x;

  // mfem::HyprePCG solver(A);
  // solver.SetTol(1e-12);
  // solver.SetMaxIter(500);
  // solver.SetPrintLevel(2);

  // mfem::HypreBoomerAMG precond;
  // precond.SetPrintLevel(0);
  // solver.SetPreconditioner(precond);

  // solver.Mult(b, x);  // solve Ax = B

  ////////////////////////////////////////////////////////
  /// Recover thermodynamic contribution
  ////////////////////////////////////////////////////////
}

/**
 * @brief Second order finite difference approximation
 *
 *
 * @tparam T
 * @return mfem::SparseMatrix
 */
template <typename T>
mfem::SparseMatrix CalphadBase<T>::get_A4linearKKS(
    const std::vector<std::tuple<std::string, std::string>> &chemicalsystem,
    const std::string &phase, const int node) {
  mfem::SparseMatrix AA;
  // TODO(cci) remove last_component
  //  Diagonal
  // elem_id corresponds to index of variable with delta_x, index of derivative
  int elem_id = 0;
  for (const auto &[elem, unit] : chemicalsystem) {
    double d2gd2x =
        this->chemical_potentials_right_x_[std::make_tuple(elem_id, node, elem, phase)] -
        this->chemical_potentials_right_x_[std::make_tuple(elem_id, node,
                                                           this->element_removed_from_ic_, phase)] +
        this->chemical_potentials_left_x_[std::make_tuple(elem_id, node, elem, phase)] -
        this->chemical_potentials_left_x_[std::make_tuple(elem_id, node,
                                                          this->element_removed_from_ic_, phase)];
    d2gd2x /= 2. * this->KKS_composition_increment_;
    AA.Add(elem_id, elem_id, d2gd2x);
    elem_id++;
  }
  // Off-diagonal
  elem_id = 0;
  int elem_jd = 1;
  for (const auto &[ielem, unit] : chemicalsystem) {
    elem_jd = 1;
    for (const auto &[jelem, unit] : chemicalsystem) {
      // d2dxdy = 0.5 d(mu_x - mu-n)/dy + 0.5 d(mu_y - mu-n)/dx
      if (ielem != jelem) {
        double d2gdxdy =
            0.5 *
                (this->chemical_potentials_right_x_[std::make_tuple(elem_id, node, jelem, phase)] -
                 this->chemical_potentials_right_x_[std::make_tuple(
                     elem_id, node, this->element_removed_from_ic_, phase)]) +
            0.5 * (this->chemical_potentials_left_x_[std::make_tuple(elem_jd, node, ielem, phase)] -
                   this->chemical_potentials_left_x_[std::make_tuple(
                       elem_jd, node, this->element_removed_from_ic_, phase)]);
        d2gdxdy /= 2. * this->KKS_composition_increment_;

        AA.Add(elem_id, elem_jd, d2gdxdy);

        elem_jd++;
      }
    }
    elem_id++;
  }
  return AA;
}

/**
 * @brief First order finite difference approximation
 *
 * @tparam T
 * @param phase
 * @param node
 * @return mfem::Vector
 */
template <typename T>
mfem::Vector CalphadBase<T>::get_m4linearKKS(
    const std::vector<std::tuple<std::string, std::string>> &chemicalsystem,
    const std::string &phase, const int node) {
  mfem::Vector mm;
  int elem_id = 0;
  for (const auto &[ielem, unit] : chemicalsystem) {
    mm(elem_id) = this->chemical_potentials_by_phase_[std::make_tuple(node, ielem, phase)] -
                  this->chemical_potentials_by_phase_[std::make_tuple(
                      node, this->element_removed_from_ic_, phase)];
    elem_id++;
  }

  return mm;
}

/**
 * @brief Second order finite difference approximation (cross derivatives)
 *
 * @tparam T
 * @param phase
 * @param node
 * @return mfem::Vector
 */
template <typename T>
mfem::Vector CalphadBase<T>::get_h4linearKKS(
    const std::vector<std::tuple<std::string, std::string>> &chemicalsystem,
    const std::string &phase, const int node) {
  mfem::Vector hh;
  int elem_id = 0;
  for (const auto &[ielem, unit] : chemicalsystem) {
    hh(elem_id) = this->chemical_potentials_right_T_[std::make_tuple(node, ielem, phase)] -
                  this->chemical_potentials_right_T_[std::make_tuple(
                      node, this->element_removed_from_ic_, phase)] +
                  this->chemical_potentials_left_T_[std::make_tuple(node, ielem, phase)] -
                  this->chemical_potentials_left_T_[std::make_tuple(
                      node, this->element_removed_from_ic_, phase)];
    hh(elem_id) /= 2. * this->KKS_temperature_increment_;
    elem_id++;
  }

  return hh;
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

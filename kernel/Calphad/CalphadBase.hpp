
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
#include "Coefficients/PhaseFieldPotentials.hpp"
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
  double KKS_nucleation_seed_;
  double KKS_threshold_;
  bool KKS_temperature_scheme_;
  bool KKS_check_nucleation_;

  PotentialFunctions<0, ThermodynamicsPotentialDiscretization::Implicit,
                     ThermodynamicsPotentials::H>
      interpolation_func_;
  void get_KKS_parameters();
  void KKS_execute(const int dt, const std::vector<T> &tp_gf, const std::vector<T> &tp_gf_old,
                   const std::tuple<std::string, T, T> &phasefields_gf,
                   const std::vector<std::tuple<std::string, std::string>> &chemicalsystem,
                   const std::vector<std::tuple<std::string, std::string, T, T>> &x_gf);

 protected:
  std::shared_ptr<CalphadUtils<T>> CU_;
  std::string description_{"UNKNOWN CALPHAD"};
  std::string element_removed_from_ic_;
  const Parameters &params_;

  // KKS
  mfem::SparseMatrix *get_A4linearKKS(
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
  std::map<std::tuple<int, std::string>, double> diffusion_chemical_potentials_;
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
      const int dt, const size_t nb_nodes,
      std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>> &output_system);

  virtual void execute(const int dt, const std::set<int> &list_nodes, const std::vector<T> &tp_gf,
                       const std::vector<std::tuple<std::string, std::string>> &chemicalsystem,
                       std::optional<std::vector<std::tuple<std::string, std::string>>>
                           status_phase = std::nullopt) = 0;

 public:
  explicit CalphadBase(const Parameters &params);
  CalphadBase(const Parameters &params, bool is_KKS);

  virtual void initialize(
      const std::vector<std::tuple<std::string, std::string>> &sorted_chemical_system) = 0;

  void global_execute(
      const int dt, const std::vector<T> &tp_gf,
      const std::vector<std::tuple<std::string, std::string>> &chemicalsystem,
      std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>> &output_system,
      std::optional<const std::tuple<std::string, T, T>> phase_fields = std::nullopt,
      std::optional<const std::vector<T>> tp_gf_old = std::nullopt,
      std::optional<const std::vector<std::tuple<std::string, std::string, T, T>>> x_gf =
          std::nullopt);

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
CalphadBase<T>::CalphadBase(const Parameters &params) : CalphadBase(params, false) {
  this->CU_ = std::make_shared<CalphadUtils<T>>();
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
CalphadBase<T>::CalphadBase(const Parameters &params, bool is_KKS)
    : params_(params), is_KKS_(is_KKS) {
  this->CU_ = std::make_shared<CalphadUtils<T>>();

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
  // this->KKS_check_nucleation_ =
  //     this->params_.template get_param_value_or_default<bool>("KKS_check_nucleation", true);
  // this->KKS_nucleation_seed_ =
  //     this->params_.template get_param_value_or_default<double>("KKS_nucleation_seed", true);

  this->KKS_secondary_phase_ =
      this->params_.template get_param_value<std::string>("KKS_secondary_phase");

  this->KKS_temperature_increment_ =
      this->params_.template get_param_value<double>("KKS_temperature_increment");
  this->KKS_composition_increment_ =
      this->params_.template get_param_value<double>("KKS_composition_increment");
  this->KKS_temperature_scheme_ = this->params_.template get_param_value_or_default<bool>(
      "KKS_temperature_explicit_scheme", false);
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
    const int dt, const std::vector<T> &tp_gf, const std::vector<T> &tp_gf_old,
    const std::tuple<std::string, T, T> &phasefields_gf,
    const std::vector<std::tuple<std::string, std::string>> &chemicalsystem,
    const std::vector<std::tuple<std::string, std::string, T, T>> &x_gf) {
  // TODO(cci): gerer le cas xalpha_old < xmin_tol ou premier dt, prendre x

  const double xmin_tol = 1.e-12;
  const int nb_elem = chemicalsystem.size();
  // Creation initial list of nodes
  const size_t nb_nodes = this->CU_->get_size(tp_gf[0]);
  std::set<int> list_nodes;
  for (int i = 0; i < nb_nodes; ++i) {
    list_nodes.insert(i);
  }
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
  std::set<int> indices_ph_newly_inter;
  for (int i = 0; i < nb_nodes; ++i) {
    if (phi_gf[i] > 1 - this->KKS_threshold_) {
      indices_ph_1.insert(i);
    } else if (phi_gf[i] < this->KKS_threshold_) {
      SlothInfo::print(" LinearKKS verification :liuqid e node  ", i);
      indices_ph_2.insert(i);
    } else {
      SlothInfo::print(" LinearKKS verification : interfacial node  ", i, " phi_gf_old[i] ",
                       phi_gf_old[i], " phi_gf[i] ", phi_gf[i]);
      indices_inter.insert(i);
      // Node previously solid -> no liquid
      if (phi_gf_old[i] > 1 - this->KKS_threshold_) {
        indices_ph_newly_inter.insert(i);
      }
    }

    // Interpolation function at node
    FType H = this->interpolation_func_.getPotentialFunction(phi_gf_old(i));
    Hphi(i) = H(phi_gf(i));
    Hphi_old(i) = H(phi_gf_old(i));
  }

  ////////////////////////////////////////////////
  // Point of expansion (Tbar,Xbar)         //////
  ////////////////////////////////////////////////
  // Phase : index _1
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

  // Xbar

  MFEM_VERIFY(x_gf.size() > 0, "Error while getting molar fraction by phase");
  for (int j = 0; j < x_gf.size(); j++) {
    const auto &[elem, phase_elem, x_elem, x_elem_old] = x_gf[j];
    // TODO(cci) improve with findIndexOfTuple
    auto it = std::find_if(
        chemicalsystem.begin(), chemicalsystem.end(),
        [&elem](const std::tuple<std::string, std::string> &t) { return std::get<0>(t) == elem; });
    int id = -1;
    if (it != chemicalsystem.end()) {
      id = std::distance(chemicalsystem.begin(), it);
    }

    MFEM_VERIFY(id > -1, "Error while getting the element");

    if (phase_elem == phase) {
      if (elem != this->element_removed_from_ic_) {
        bar_tp_gf_ph_1[id + 2] = x_elem;
      } else {
        bar_tp_gf_ph_1[id + 2] = tp_gf[id + 2];
      }
      pure_bar_tp_gf_ph_1[id + 2] = tp_gf[id + 2];

    } else if (phase_elem == this->KKS_secondary_phase_) {
      if (elem != this->element_removed_from_ic_) {
        bar_tp_gf_ph_2[id + 2] = x_elem;

        // When liquid is previously absent, initialize with solid
        if (!indices_inter.empty()) {
          for (const int id_new : indices_inter) {
            if (bar_tp_gf_ph_2[id + 2](id_new) < xmin_tol) {
              std::cout << " bar_tp_gf_ph_2 : new " << id_new << std::endl;
              bar_tp_gf_ph_2[id + 2](id_new) = tp_gf[id + 2](id_new);
            }
          }
        }
      } else {
        bar_tp_gf_ph_2[id + 2] = tp_gf[id + 2];
      }
      pure_bar_tp_gf_ph_2[id + 2] = tp_gf[id + 2];

    } else {
      std::runtime_error(
          "Error while setting molar fraction at point approximation. Unknown phase " + phase_elem);
    }
  }
  for (const auto &vv : bar_tp_gf_ph_1) {
    std::cout << " bar_tp_gf_ph_1 " << vv.Size() << " " << vv(0) << std::endl;
  }
  for (const auto &vv : bar_tp_gf_ph_2) {
    std::cout << " bar_tp_gf_ph_2 " << vv.Size() << " " << vv(0) << std::endl;
  }

  ////////////////////////////////////////////////////////
  /// Calculation of all thermodynamic contributions
  ////////////////////////////////////////////////////////

  // Interface calculations
  auto calculate_interface =
      [&](const std::set<int> &indices, const std::vector<T> &delta_tp_gf, const double increment,
          const int id_incr,
          std::map<std::tuple<int, std::string, std::string>, double> &chemical_potential_interface,
          std::vector<std::tuple<std::string, std::string>> given_phase,
          const std::string &primary_phase) {
        std::vector<T> delta_tp = delta_tp_gf;
        if (id_incr > -1) {
          delta_tp[id_incr] += increment;
        }
        this->execute(dt, indices, delta_tp, chemicalsystem, given_phase);

        for (const auto &in : indices) {
          for (const auto &elem : chemicalsystem) {
            const auto &[elem1, unit] = elem;
            const double mu = this->chemical_potentials_.at(std::make_tuple(in, elem1));
            chemical_potential_interface.emplace(std::make_tuple(in, elem1, primary_phase), mu);
          }
        }
      };
  //
  auto calculate_interface_x = [&](const std::vector<T> &delta_tp_gf, const double increment,
                                   const int id_incr, int index_el,
                                   std::map<std::tuple<int, int, std::string, std::string>, double>
                                       &chemical_potential_interface,
                                   std::vector<std::tuple<std::string, std::string>> given_phase,
                                   const std::string &primary_phase) {
    std::vector<T> delta_tp = delta_tp_gf;
    if (id_incr > -1) {
      delta_tp[id_incr] += increment;
    }

    this->execute(dt, indices_inter, delta_tp, chemicalsystem, given_phase);
    for (const auto &in : indices_inter) {
      for (const auto &elem : chemicalsystem) {
        const auto &[elem1, unit] = elem;
        const double mu = this->chemical_potentials_.at(std::make_tuple(in, elem1));
        chemical_potential_interface.emplace(std::make_tuple(index_el, in, elem1, primary_phase),
                                             mu);
      }
    }
  };
  // TODO(cci) reordonner pour ne pas faire inutilement sur les phases non presentes et sur
  // interface
  // TODO(cci) : nucleation faire un test sur l'apparition de la secondary_phase
  //  T : equilibrium calculations in both phases (Tbar, Xbar)
  // calculate_interface(indices_ph_1, bar_tp_gf_ph_1, phase,
  // this->chemical_potentials_by_phase_); For primary phase it is assumed that another phase can
  // appear
  std::vector<std::tuple<std::string, std::string>> st_phase_12 = {
      {phase, "entered"}, {this->KKS_secondary_phase_, "dormant"}};
  std::vector<std::tuple<std::string, std::string>> st_phase_1 = {{phase, "entered"}};
  std::vector<std::tuple<std::string, std::string>> st_phase_2 = {
      {this->KKS_secondary_phase_, "entered"}};

  calculate_interface(indices_ph_1, pure_bar_tp_gf_ph_1, 0., -1,
                      this->chemical_potentials_by_phase_, st_phase_12, phase);

  // Check "Nucleation" with driving force
  for (const auto &[key, dgm] : this->driving_forces_) {
    const auto &[node, phase] = key;
    if (dgm > 0.) {
      std::cout << " Nucleation start at node " << node << " for phase " << phase << std::endl;
    } else {
      // ???? TODO(cci) := to be confirmed
      this->driving_forces_[std::make_tuple(node, phase)] = 0.;
    }
  }
  ///////
  // If node of secondary phase exist
  if (!indices_ph_2.empty()) {
    SlothInfo::print(" LinearKKS verification : check calculation in LIQUID  ");

    calculate_interface(indices_ph_2, pure_bar_tp_gf_ph_2, 0., -1,
                        this->chemical_potentials_by_phase_, st_phase_2,
                        this->KKS_secondary_phase_);
    SlothInfo::print(" LinearKKS verification : check calculation in LIQUID end  ");
  }
  // If node of interface exist
  if (!indices_inter.empty()) {
    for (const auto &vv : bar_tp_gf_ph_1) {
      std::cout << "solid inyterface  bar_tp_gf_ph_1 " << vv.Size() << " " << vv(0) << std::endl;
    }
    for (const auto &vv : bar_tp_gf_ph_2) {
      std::cout << "liquid interface  bar_tp_gf_ph_2 " << vv.Size() << " " << vv(0) << std::endl;
    }

    // T-dT : equilibrium calculations in both phases (Tbar - deltaT, Xbar)
    SlothInfo::print(" LinearKKS verification : check calculation in interface 0 ");

    calculate_interface(indices_inter, bar_tp_gf_ph_1, -this->KKS_temperature_increment_, 0,
                        this->chemical_potentials_left_T_, st_phase_1, phase);

    SlothInfo::print(" LinearKKS verification : check calculation in interface 1 ");
    calculate_interface(indices_inter, bar_tp_gf_ph_2, -this->KKS_temperature_increment_, 0,
                        this->chemical_potentials_left_T_, st_phase_2, this->KKS_secondary_phase_);

    // T+dT : equilibrium calculations in both phases (Tbar + deltaT, Xbar)

    SlothInfo::print(" LinearKKS verification : check calculation in interface 2 ");
    calculate_interface(indices_inter, bar_tp_gf_ph_1, this->KKS_temperature_increment_, 0,
                        this->chemical_potentials_right_T_, st_phase_1, phase);

    SlothInfo::print(" LinearKKS verification : check calculation in interface 3 ");
    calculate_interface(indices_inter, bar_tp_gf_ph_2, this->KKS_temperature_increment_, 0,
                        this->chemical_potentials_right_T_, st_phase_2, this->KKS_secondary_phase_);

    SlothInfo::print(" LinearKKS verification : check calculation in interface 4 ");
    // Loop over chemical system (except the reference element)
    int ielem = 0;
    for (const auto &[elem, unit] : chemicalsystem) {
      if (elem != this->element_removed_from_ic_) {
        // x + dx
        // Equilibrium calculations in both phases at (Tbar, Xbar + delta X)

        SlothInfo::print(" LinearKKS verification : check calculation in interface 5 ", elem);
        calculate_interface_x(bar_tp_gf_ph_1, this->KKS_composition_increment_, ielem + 2, ielem,
                              this->chemical_potentials_right_x_, st_phase_1, phase);

        SlothInfo::print(" LinearKKS verification : check calculation in interface 6 ", elem);
        calculate_interface_x(bar_tp_gf_ph_2, this->KKS_composition_increment_, ielem + 2, ielem,
                              this->chemical_potentials_right_x_, st_phase_2,
                              this->KKS_secondary_phase_);

        // x - dx
        // Equilibrium calculations in both phases at (Tbar, Xbar - delta X)

        SlothInfo::print(" LinearKKS verification : check calculation in interface 7 ", elem);
        calculate_interface_x(bar_tp_gf_ph_1, -this->KKS_composition_increment_, ielem + 2, ielem,
                              this->chemical_potentials_left_x_, st_phase_1, phase);

        SlothInfo::print(" LinearKKS verification : check calculation in interface 8 ", elem);
        calculate_interface_x(bar_tp_gf_ph_2, -this->KKS_composition_increment_, ielem + 2, ielem,
                              this->chemical_potentials_left_x_, st_phase_2,
                              this->KKS_secondary_phase_);
      }
      ielem++;
    }

    ////////////////////////////////////////////////////////
    /// Solve linear systems
    ////////////////////////////////////////////////////////
    SlothInfo::print(" LinearKKS verification : check 0  ");
    for (const auto &node : indices_inter) {
      SlothInfo::print(" LinearKKS verification : check 1  ");
      mfem::SparseMatrix *Al =
          this->get_A4linearKKS(chemicalsystem, this->KKS_secondary_phase_, node);
      SlothInfo::print(" LinearKKS verification : check 1 AL ");
      SlothInfo::print(" LinearKKS verification : check 1.1  ");
      mfem::SparseMatrix *As = this->get_A4linearKKS(chemicalsystem, phase, node);

      SlothInfo::print(" LinearKKS verification : check 1.1 AS ");
      mfem::Vector hs = this->get_h4linearKKS(chemicalsystem, phase, node);
      SlothInfo::print(" LinearKKS verification : check 1.2 hs ");
      mfem::Vector hl = this->get_h4linearKKS(chemicalsystem, this->KKS_secondary_phase_, node);
      SlothInfo::print(" LinearKKS verification : check 1.3 hl ");
      mfem::Vector ms = this->get_m4linearKKS(chemicalsystem, phase, node);
      SlothInfo::print(" LinearKKS verification : check 1.4 ms ");
      mfem::Vector ml = this->get_m4linearKKS(chemicalsystem, this->KKS_secondary_phase_, node);
      SlothInfo::print(" LinearKKS verification : check 1.5 ml ");
      SlothInfo::print(" LinearKKS verification : check 2  ");
      //================
      // ml - ms + (T - Tbar) * (hl - hs)
      // ml - ms
      mfem::Vector ml_minus_ms = ml;
      ml_minus_ms -= ms;
      // T - Tbar
      mfem::Vector deltaT = tp_gf[0];
      deltaT -= bar_tp_gf_ph_1[0];
      // hl * (T - Tbar)
      mfem::Vector hl_deltaT = hl;
      hl_deltaT *= deltaT;
      // hs * (T - Tbar)
      mfem::Vector hs_deltaT = hs;
      hs_deltaT *= deltaT;
      // hl * (T - Tbar) - hs * (T - Tbar)
      mfem::Vector hl_minus_hs = hl_deltaT;
      hl_minus_hs -= hs_deltaT;

      SlothInfo::print(" LinearKKS verification : check 3  ");
      // DeltaX
      mfem::Vector deltaX(nb_elem - 1);
      int ielem = 0;
      int jelem = 0;
      for (const auto &[elem, unit] : chemicalsystem) {
        if (elem != this->element_removed_from_ic_) {
          const double dx = tp_gf[jelem + 2](node) - tp_gf_old[jelem + 2](node);
          const double xdH = (bar_tp_gf_ph_1[jelem + 2](node) - bar_tp_gf_ph_2[jelem + 2](node)) *
                             (Hphi(node) - Hphi_old(node));

          deltaX(ielem) = dx - xdH;
          ielem++;
        }
        jelem++;
      }

      SlothInfo::print(" LinearKKS verification : check 4  ");
      // Assemble

      mfem::Array<int> offsets(2);
      offsets[0] = 0;
      offsets[1] = nb_elem;

      mfem::BlockVector bb(offsets);
      mfem::Vector &b0 = bb.GetBlock(0);
      b0 = ml_minus_ms;
      b0 += hl_minus_hs;
      mfem::Vector &b1 = bb.GetBlock(1);
      b1 = deltaX;

      mfem::BlockMatrix *AA = new mfem::BlockMatrix(offsets, offsets);
      AA->SetBlock(0, 0, As);
      *Al *= -1.;
      AA->SetBlock(0, 0, Al);

      // identity
      // TODO check nb_elem-1
      mfem::Vector HphiNode(nb_elem - 1);
      HphiNode = Hphi(node);
      mfem::Vector one_minus_HphiNode(nb_elem - 1);
      one_minus_HphiNode = 1. - Hphi(node);

      mfem::SparseMatrix *A10 = new mfem::SparseMatrix(HphiNode);
      mfem::SparseMatrix *A11 = new mfem::SparseMatrix(one_minus_HphiNode);
      AA->SetBlock(1, 0, A10);
      AA->SetBlock(1, 1, A11);

      AA->Finalize();
      AA->Print();
      SlothInfo::print(" LinearKKS verification : check 5  ");

      // Solve
      auto solver = std::make_shared<mfem::CGSolver>();
      auto prec = std::make_shared<mfem::BlockDiagonalPreconditioner>(offsets);
      auto smoother0 = std::make_shared<mfem::GSSmoother>(AA->GetBlock(0, 0));
      auto smoother1 = std::make_shared<mfem::GSSmoother>(AA->GetBlock(1, 1));
      prec->SetDiagonalBlock(0, smoother0.get());
      prec->SetDiagonalBlock(1, smoother1.get());

      solver->SetOperator(*AA);
      solver->SetAbsTol(1e-12);
      solver->SetRelTol(1e-12);
      solver->SetMaxIter(50);
      solver->SetPrintLevel(2);
      solver->SetPreconditioner(*prec);

      mfem::BlockVector deltaX_phase(offsets);
      solver->Mult(bb, deltaX_phase);
      // Back to Al
      *Al *= -1.;
      ////////////////////////////////////////////////////////
      /// Recover thermodynamic contribution
      ////////////////////////////////////////////////////////
      const mfem::Vector delta_XS = deltaX_phase.GetBlock(0);
      const mfem::Vector delta_XL = deltaX_phase.GetBlock(1);
      // Xs, Xl
      int i = 0;
      double sumS = 0.;
      double sumL = 0.;
      for (const auto &[elem, unit] : chemicalsystem) {
        if (elem != this->element_removed_from_ic_) {
          this->elem_mole_fraction_by_phase_[std::make_tuple(node, phase, elem)] =
              bar_tp_gf_ph_1[i + 2](node) + delta_XS(i);
          this->elem_mole_fraction_by_phase_[std::make_tuple(
              node, this->KKS_secondary_phase_, elem)] = bar_tp_gf_ph_2[i + 2](node) + delta_XL(i);
          sumS += this->elem_mole_fraction_by_phase_[std::make_tuple(node, phase, elem)];
          sumL += this->elem_mole_fraction_by_phase_[std::make_tuple(
              node, this->KKS_secondary_phase_, elem)];
        }
        ++i;
      }
      // TODO(cci) initialization with XX global to have continuity?
      this->elem_mole_fraction_by_phase_[std::make_tuple(
          node, phase, this->element_removed_from_ic_)] = 1. - sumS;
      this->elem_mole_fraction_by_phase_[std::make_tuple(
          node, this->KKS_secondary_phase_, this->element_removed_from_ic_)] = 1. - sumL;

      // Difference of chemical potentials
      mfem::Vector dmu_s = ms;
      dmu_s += hs_deltaT;
      mfem::Vector AsDeltaXs;
      As->Mult(delta_XS, AsDeltaXs);
      dmu_s += AsDeltaXs;
      mfem::Vector dmu_l = ml;
      dmu_l += hl_deltaT;
      mfem::Vector AlDeltaXl;
      Al->Mult(delta_XL, AlDeltaXl);
      dmu_l += AlDeltaXl;

      // Gibbs energies (implicit case)
      double gs = this->energies_of_phases_[std::make_tuple(node, phase, "GM")];
      double gl =
          this->energies_of_phases_[std::make_tuple(node, this->KKS_secondary_phase_, "GM")];
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
      this->energies_of_phases_[std::make_tuple(node, phase, "GM")] = gs;
      this->energies_of_phases_[std::make_tuple(node, this->KKS_secondary_phase_, "GM")] = gl;

      // Driving force by phase : contribution -gm
      this->driving_forces_[std::make_tuple(node, phase)] =
          -this->energies_of_phases_[std::make_tuple(node, phase, "GM")];
      this->driving_forces_[std::make_tuple(node, this->KKS_secondary_phase_)] =
          -this->energies_of_phases_[std::make_tuple(node, this->KKS_secondary_phase_, "GM")];

      i = 0;
      for (const auto &[elem, unit] : chemicalsystem) {
        if (elem == this->element_removed_from_ic_) continue;
        SlothInfo::print(" LinearKKS verification at node  ", node, " elem ", elem,
                         " DeltaX to obtain ", tp_gf[i + 2](node) - tp_gf_old[i + 2](node),
                         " DeltaX predicted ",
                         delta_XS(i) * Hphi(node) + delta_XL(i) * (1. - Hphi(node)));
        SlothInfo::print(
            " LinearKKS verification at node  ", node, "  elem ", elem, " X to obtain ",
            tp_gf[i + 2](node), " X predicted ",
            this->elem_mole_fraction_by_phase_[std::make_tuple(node, phase, elem)] * Hphi(node) +
                this->elem_mole_fraction_by_phase_[std::make_tuple(node, this->KKS_secondary_phase_,
                                                                   elem)] *
                    (1. - Hphi(node)));

        SlothInfo::print(" LinearKKS verification at node : dmu  ", node, " elem ", i, " dmu_s ",
                         dmu_s(i), " dmu_l ", dmu_l(i));
        // Diffusion chemical potential
        this->diffusion_chemical_potentials_[std::make_tuple(node, elem)] = dmu_s(i);

        // Driving force : contribution mu_i * x_i
        this->driving_forces_[std::make_tuple(node, phase)] +=
            this->diffusion_chemical_potentials_[std::make_tuple(node, elem)] *
            this->elem_mole_fraction_by_phase_[std::make_tuple(node, phase, elem)];

        this->driving_forces_[std::make_tuple(node, this->KKS_secondary_phase_)] +=
            this->diffusion_chemical_potentials_[std::make_tuple(node, elem)] *
            this->elem_mole_fraction_by_phase_[std::make_tuple(node, this->KKS_secondary_phase_,
                                                               elem)];

        i++;
      }

      SlothInfo::print(" LinearKKS verification at node g: ", node, " gs ", gs, " gl ", gl);

      SlothInfo::print(" LinearKKS verification at node dgm: ", node, " dgm_s ",
                       this->driving_forces_[std::make_tuple(node, phase)], " dgm_l ",
                       this->driving_forces_[std::make_tuple(node, this->KKS_secondary_phase_)]);
    }
  }
}

/**
 * @brief Second order finite difference approximation
 *
 *
 * @tparam T
 * @return mfem::SparseMatrix
 */
template <typename T>
mfem::SparseMatrix *CalphadBase<T>::get_A4linearKKS(
    const std::vector<std::tuple<std::string, std::string>> &chemicalsystem,
    const std::string &phase, const int node) {
  const int nb_elem = chemicalsystem.size();
  mfem::SparseMatrix *AA = new mfem::SparseMatrix(nb_elem - 1, nb_elem - 1);

  SlothInfo::print(" get_A4linearKKS ");
  // TODO(cci) remove last_component
  //  Diagonal
  // elem_id corresponds to index of variable with delta_x, index of derivative
  int elem_id = 0;
  for (const auto &[elem, unit] : chemicalsystem) {
    if (elem != this->element_removed_from_ic_) {
      double d2gd2x =
          this->chemical_potentials_right_x_[std::make_tuple(elem_id, node, elem, phase)] -
          this->chemical_potentials_right_x_[std::make_tuple(
              elem_id, node, this->element_removed_from_ic_, phase)] +
          this->chemical_potentials_left_x_[std::make_tuple(elem_id, node, elem, phase)] -
          this->chemical_potentials_left_x_[std::make_tuple(elem_id, node,
                                                            this->element_removed_from_ic_, phase)];
      d2gd2x /= 2. * this->KKS_composition_increment_;
      std::cout << " diag elem_id " << elem_id << d2gd2x << " elem " << elem << std::endl;
      AA->Set(elem_id, elem_id, d2gd2x);
      std::cout << " diag elem_id end " << std::endl;
      elem_id++;
    }
  }
  // Off-diagonal
  elem_id = 0;
  int elem_jd = 1;
  for (const auto &[ielem, unit] : chemicalsystem) {
    if (ielem != this->element_removed_from_ic_) {
      elem_jd = 1;
      for (const auto &[jelem, unit] : chemicalsystem) {
        if (jelem != this->element_removed_from_ic_) {
          // d2dxdy = 0.5 d(mu_x - mu-n)/dy + 0.5 d(mu_y - mu-n)/dx
          if (ielem != jelem) {
            double d2gdxdy = 0.5 * (this->chemical_potentials_right_x_[std::make_tuple(
                                        elem_id, node, jelem, phase)] -
                                    this->chemical_potentials_right_x_[std::make_tuple(
                                        elem_id, node, this->element_removed_from_ic_, phase)]) +
                             0.5 * (this->chemical_potentials_left_x_[std::make_tuple(
                                        elem_jd, node, ielem, phase)] -
                                    this->chemical_potentials_left_x_[std::make_tuple(
                                        elem_jd, node, this->element_removed_from_ic_, phase)]);
            d2gdxdy /= 2. * this->KKS_composition_increment_;
            std::cout << " off diag  elem_id " << elem_id << " elem_jd " << elem_jd << " d2gdxdy  "
                      << d2gdxdy << " elem " << ielem << " jelem " << jelem << std::endl;

            AA->Set(elem_id, elem_jd, d2gdxdy);
            std::cout << " off diag  end " << std::endl;

            elem_jd++;
          }
        }
      }
      elem_id++;
    }
  }

  SlothInfo::print(" get_A4linearKKS 1 ");

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
  const int nb_elem = chemicalsystem.size();
  mfem::Vector mm(nb_elem - 1);
  int elem_id = 0;
  for (const auto &[ielem, unit] : chemicalsystem) {
    if (ielem != this->element_removed_from_ic_) {
      mm(elem_id) = this->chemical_potentials_by_phase_[std::make_tuple(node, ielem, phase)] -
                    this->chemical_potentials_by_phase_[std::make_tuple(
                        node, this->element_removed_from_ic_, phase)];
      elem_id++;
    }
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
  const int nb_elem = chemicalsystem.size();
  mfem::Vector hh(nb_elem - 1);
  int elem_id = 0;
  for (const auto &[ielem, unit] : chemicalsystem) {
    if (ielem != this->element_removed_from_ic_) {
      hh(elem_id) = this->chemical_potentials_right_T_[std::make_tuple(node, ielem, phase)] -
                    this->chemical_potentials_right_T_[std::make_tuple(
                        node, this->element_removed_from_ic_, phase)] +
                    this->chemical_potentials_left_T_[std::make_tuple(node, ielem, phase)] -
                    this->chemical_potentials_left_T_[std::make_tuple(
                        node, this->element_removed_from_ic_, phase)];
      hh(elem_id) /= 2. * this->KKS_temperature_increment_;
      elem_id++;
    }
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
    std::optional<const std::tuple<std::string, T, T>> phase_field_gf,
    std::optional<const std::vector<T>> tp_gf_old,
    std::optional<const std::vector<std::tuple<std::string, std::string, T, T>>> x_gf) {
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
    MFEM_VERIFY(tp_gf_old.has_value(), "Error: tp_gf_old is required for KKS execution.");
    MFEM_VERIFY(x_gf.has_value(), "Error: x_gf is required for KKS execution.");
    this->KKS_execute(dt, tp_gf, *tp_gf_old, *phase_field_gf, chemicalsystem, *x_gf);
  }
  // Use containers to update output_system
  this->update_outputs(dt, nb_nodes, output_system);
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
  this->chemical_potentials_by_phase_.clear();
  this->chemical_potentials_left_T_.clear();
  this->chemical_potentials_left_x_.clear();
  this->chemical_potentials_right_T_.clear();
  this->chemical_potentials_right_x_.clear();
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
    const int dt, const size_t nb_nodes,
    std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>> &output_system) {
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

  for (auto &[output_infos, output_value] : output_system) {
    const std::string &output_type = output_infos.back();

    // Fill output with the relevant values
    switch (calphad_outputs::from(output_type)) {
      case calphad_outputs::mu: {
        const std::string &output_elem = output_infos[1];
        for (std::size_t i = 0; i < nb_nodes; ++i) {
          output[i] =
              get_or_default(this->chemical_potentials_, std::make_tuple(i, output_elem), 0.);
        }
        break;
      }
      case calphad_outputs::dmu: {
        const std::string &output_elem = output_infos[1];
        for (std::size_t i = 0; i < nb_nodes; ++i) {
          output[i] = get_or_default(this->diffusion_chemical_potentials_,
                                     std::make_tuple(i, output_elem), 0.);
        }
        break;
      }
      case calphad_outputs::xp: {
        const std::string &output_elem = output_infos[1];
        const std::string &output_phase = output_infos[2];
        for (std::size_t i = 0; i < nb_nodes; ++i) {
          output[i] = get_or_default(this->elem_mole_fraction_by_phase_,
                                     std::make_tuple(i, output_phase, output_elem), 0.);
        }
        break;
      }
      case calphad_outputs::y: {
        const std::string &output_cons = output_infos[1];
        const int &output_sub = std::stoi(output_infos[2]);
        const std::string &output_phase = output_infos[3];
        for (std::size_t i = 0; i < nb_nodes; ++i) {
          output[i] = get_or_default(this->site_fraction_,
                                     std::make_tuple(i, output_phase, output_cons, output_sub), 0.);
        }
        break;
      }
      case calphad_outputs::g: {
        const std::string &output_phase = output_infos[1];
        for (std::size_t i = 0; i < nb_nodes; ++i) {
          output[i] =
              get_or_default(this->energies_of_phases_, std::make_tuple(i, output_phase, "G"), 0.);
        }
        break;
      }
      case calphad_outputs::gm: {
        const std::string &output_phase = output_infos[1];
        for (std::size_t i = 0; i < nb_nodes; ++i) {
          output[i] =
              get_or_default(this->energies_of_phases_, std::make_tuple(i, output_phase, "GM"), 0.);
        }
        break;
      }
      case calphad_outputs::h: {
        const std::string &output_phase = output_infos[1];
        for (std::size_t i = 0; i < nb_nodes; ++i) {
          output[i] =
              get_or_default(this->energies_of_phases_, std::make_tuple(i, output_phase, "H"), 0.);
        }
        break;
      }
      case calphad_outputs::hm: {
        const std::string &output_phase = output_infos[1];
        for (std::size_t i = 0; i < nb_nodes; ++i) {
          output[i] =
              get_or_default(this->energies_of_phases_, std::make_tuple(i, output_phase, "HM"), 0.);
        }
        break;
      }
      case calphad_outputs::dgm: {
        const std::string &output_phase = output_infos[1];
        for (std::size_t i = 0; i < nb_nodes; ++i) {
          output[i] = get_or_default(this->driving_forces_, std::make_tuple(i, output_phase), 0.);
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
          output[i] =
              get_or_default(this->mobilities_, std::make_tuple(i, output_phase, output_elem),
                             -std::numeric_limits<double>::max());
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

/**
 * @file AnalyticalParaboloidForTwoPhase.hpp
 * @author cp273896  (clement.plumecocq@cea.fr)
 * @brief Analytical thermodynamic description for a two-phase problem
 * @version 0.1
 * @date 2025-01-15
 *
 * Copyright CEA (c) 2025
 *
 */
#include <algorithm>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "Calphad/CalphadBase.hpp"
#include "Calphad/CalphadUtils.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"

#pragma once

class NonLinearSys_CALPHAD : public mfem::Operator {
 public:
  mutable mfem::DenseMatrix J;
  NonLinearSys_CALPHAD() : Operator(2), J(2, 2) {}

  mfem::real_t varphi, c, c_eq_s, c_eq_l;
  double interpolationp(double varphi) const {
    // return std::pow(varphi,3) * (10 - 15 * varphi + 6 * std::pow(varphi,2));
    return varphi;
  }
  double mu(double c, double ceq) const { return ((c - ceq)); }

  virtual void Mult(const mfem::Vector& x, mfem::Vector& y) const override {
    double x1 = x(0);
    double x2 = x(1);

    y(0) = interpolationp(this->varphi) * x1 + (1 - interpolationp(this->varphi)) * x2 - this->c;
    y(1) = mu(x1, c_eq_s) - mu(x2, c_eq_l);
  }

  mfem::Operator& GetGradient(const mfem::Vector& x) const override {
    double x1 = x(0);
    double x2 = x(1);
    double eps = 1e-10;

    J(0, 0) = interpolationp(this->varphi);
    J(0, 1) = 1 - interpolationp(this->varphi);
    J(1, 0) = 1.;
    J(1, 1) = -1.;

    return J;
  }

  void get_infos(const double& varphi, std::vector<double>& c_, const double& c_eq_s_,
                 const double& c_eq_l_) {
    this->varphi = varphi;
    this->c = c_[0];
    this->c_eq_s = c_eq_s_;
    this->c_eq_l = c_eq_l_;
  }

};  // end class

template <typename T>
class AnalyticalParaboloidForTwoPhase : CalphadBase<T> {
 private:
  std::map<std::tuple<std::string, std::string>, mfem::real_t> k_coeff, c_eq;
  std::unique_ptr<CalphadUtils<T>> CU_;
  void compute(
      const size_t nb_nodes, const std::vector<T>& tp_gf,
      const std::vector<std::tuple<std::string, std::string>>& chemical_system,
      std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system);

  void check_variables_consistency(
      std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system);

 public:
  explicit AnalyticalParaboloidForTwoPhase(const Parameters& params);

  void initialize() override;

  void execute(const int dt, const std::vector<T>& aux_gf,
               const std::vector<std::tuple<std::string, std::string>>& chemical_system,
               std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>&
                   output_system) override;

  void finalize() override;

  ////////////////////////////////

  void get_parameters() override;

  mfem::real_t get_energy(mfem::real_t x, mfem::real_t val_eq, mfem::real_t k = 1.) {
    return 0.5 * k * std::pow((x - val_eq), 2);
  };

  mfem::real_t get_mu(mfem::real_t x, mfem::real_t val_eq, mfem::real_t k = 1.) {
    return k * (x - val_eq);
  };

  ////////////////////////////////

  ~AnalyticalParaboloidForTwoPhase();
};

////////////////////////////////
////////////////////////////////

/**
 * @brief Get the parameters associated with the AnalyticalParaboloidForTwoPhase object
 *
 * @tparam T
 */
template <typename T>
void AnalyticalParaboloidForTwoPhase<T>::get_parameters() {
  this->description_ = this->params_.template get_param_value_or_default<std::string>(
      "description", "Analytical thermodynamic description for an ideal solution ");
  this->k_coeff =
      this->params_
          .template get_param_value<std::map<std::tuple<std::string, std::string>, mfem::real_t>>(
              "coefficient_k");
  this->c_eq =
      this->params_
          .template get_param_value<std::map<std::tuple<std::string, std::string>, mfem::real_t>>(
              "equilibrium_composition");
}

////////////////////////////////
////////////////////////////////
/**
 * @brief Construct a new AnalyticalParaboloidForTwoPhase::AnalyticalParaboloidForTwoPhase object
 *
 * @param params
 */
template <typename T>
AnalyticalParaboloidForTwoPhase<T>::AnalyticalParaboloidForTwoPhase(const Parameters& params)
    : CalphadBase<T>(params) {
  this->CU_ = std::make_unique<CalphadUtils<T>>();
  this->get_parameters();
}

/**
 * @brief Initialization of the thermodynamic calculation
 *
 * @tparam T
 */
template <typename T>
void AnalyticalParaboloidForTwoPhase<T>::initialize() {}

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
void AnalyticalParaboloidForTwoPhase<T>::execute(
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
void AnalyticalParaboloidForTwoPhase<T>::compute(
    const size_t nb_nodes, const std::vector<T>& tp_gf,
    const std::vector<std::tuple<std::string, std::string>>& chemical_system,
    std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system) {
  // Let us assume an ideal mixing solution
  const std::vector<std::string>& phase = {"LIQUID", "SOLID"};
  const std::vector<std::string> energy_names = {"G"};
  std::vector<double> tp_gf_at_node(tp_gf.size());

  const auto temperature_sort_method =
      this->params_.template get_param_value_or_default<std::string>("temperature_sort_method",
                                                                     "No");
  const auto pressure_sort_method =
      this->params_.template get_param_value_or_default<std::string>("pressure_sort_method", "No");

  std::vector<int> sorted_n_t_p =
      this->CU_->sort_nodes(tp_gf[0], tp_gf[1], temperature_sort_method, pressure_sort_method);

  // Process CALPHAD calculations for each node
  for (const auto& id : sorted_n_t_p) {
    // Populate tp_gf_at_node for the current node
    std::transform(tp_gf.begin(), tp_gf.end(), tp_gf_at_node.begin(),
                   [&id](const T& vec) { return vec[id]; });
    const auto Temp = tp_gf_at_node[0];
    const auto x = tp_gf_at_node[3];
    const auto varphi = tp_gf_at_node[2];
    const int index_ref_x = 3;

    bool linear_solve = true;
    std::vector<std::tuple<std::string, double>> tuple_x;
    for (size_t i = index_ref_x; i < tp_gf_at_node.size(); ++i) {
      tuple_x.push_back(
          std::make_tuple(std::get<0>(chemical_system[i - index_ref_x]), tp_gf_at_node[i]));
    }
    this->energies_of_phases_[std::make_tuple(id, "SOLID", energy_names[0])] = 0.;
    this->energies_of_phases_[std::make_tuple(id, "LIQUID", energy_names[0])] = 0.;

    for (const auto& element : tuple_x) {
      std::string current_elem = std::get<0>(element);
      mfem::real_t current_x = std::get<1>(element);
      mfem::real_t k_solid = this->k_coeff[std::make_tuple("SOLID", current_elem)];
      mfem::real_t k_liquid = this->k_coeff[std::make_tuple("LIQUID", current_elem)];
      mfem::real_t c_eq_l = this->c_eq[std::make_tuple("LIQUID", current_elem)];
      mfem::real_t c_eq_s = this->c_eq[std::make_tuple("SOLID", current_elem)];
      this->chemical_potentials_[std::make_tuple(id, current_elem)] =
          (varphi < 0.5) ? this->get_mu(current_x, c_eq_l, k_liquid)
                         : this->get_mu(current_x, c_eq_s, k_solid);
      for (const auto& p : phase) {
        // Molar fraction
        const auto& key_mole_fraction_by_phase = std::make_tuple(id, p, current_elem);
        double crit = 1e-5;
        if (varphi < crit || varphi > (1 - crit)) {
          this->elem_mole_fraction_by_phase_[key_mole_fraction_by_phase] = current_x;
          if (varphi > 0.5) {
            this->elem_mole_fraction_by_phase_[std::make_tuple(id, "SOLID", current_elem)] =
                current_x;
            this->elem_mole_fraction_by_phase_[std::make_tuple(id, "LIQUID", current_elem)] = 0.;
          }else {
            this->elem_mole_fraction_by_phase_[std::make_tuple(id, "SOLID", current_elem)] = 0.;
            this->elem_mole_fraction_by_phase_[std::make_tuple(id, "LIQUID", current_elem)] =
                current_x;
          }  

        } else {
          // this->interface_resolution(this->elem_mole_fraction_by_phase_[key_mole_fraction_by_phase],varphi);
          if (linear_solve) {
            mfem::real_t p_phi = varphi;
            mfem::real_t deno = 1 - (p_phi - 1) / p_phi; 
            mfem::real_t temp1 = current_x / p_phi + (p_phi-1) / p_phi * (c_eq_l - c_eq_s);
            this->elem_mole_fraction_by_phase_[std::make_tuple(id, "SOLID", current_elem)] =  temp1 / deno;
            this->elem_mole_fraction_by_phase_[std::make_tuple(id, "LIQUID", current_elem)] = this->elem_mole_fraction_by_phase_[std::make_tuple(id, "SOLID", current_elem)] - c_eq_s + c_eq_l;
            // if (varphi < 0.5) {
            //   // solve c_s = (p-1)/p c_l + c / p
            //   mfem::real_t p_phi = varphi;
            //   this->elem_mole_fraction_by_phase_[std::make_tuple(id, "LIQUID", current_elem)] =
            //       (current_x * k_solid + c_eq_l * k_liquid * p_phi - c_eq_s * k_solid * p_phi) /
            //       (k_liquid * p_phi - k_solid * p_phi + k_solid);
            //   this->elem_mole_fraction_by_phase_[std::make_tuple(id, "SOLID", current_elem)] =
            //       (current_x * k_liquid + c_eq_l * k_liquid * p_phi - c_eq_l * k_liquid -
            //        c_eq_s * k_solid * p_phi + c_eq_s * k_solid) /
            //       (k_liquid * p_phi - k_solid * p_phi + k_solid);
            //   mfem::real_t res = p_phi * this->elem_mole_fraction_by_phase_[std::make_tuple(
            //                                  id, "SOLID", current_elem)] +
            //                      (1 - p_phi) * this->elem_mole_fraction_by_phase_[std::make_tuple(
            //                                        id, "LIQUID", current_elem)] -
            //                      current_x;
            //   ;
            // } else {
            //   mfem::real_t p_phi = varphi;
            //   this->elem_mole_fraction_by_phase_[std::make_tuple(id, "LIQUID", current_elem)] =
            //       (c_eq_l * k_liquid * p_phi - c_eq_s * k_solid * p_phi + current_x * k_solid) /
            //       (k_liquid * p_phi - k_solid * p_phi + k_solid);
            //   this->elem_mole_fraction_by_phase_[std::make_tuple(id, "SOLID", current_elem)] =
            //       (c_eq_l * k_liquid * p_phi - c_eq_l * k_liquid - c_eq_s * k_solid * p_phi +
            //        c_eq_s * k_solid + current_x * k_liquid) /
            //       (k_liquid * p_phi - k_solid * p_phi + k_solid);
            //   mfem::real_t res = p_phi * this->elem_mole_fraction_by_phase_[std::make_tuple(
            //                                  id, "SOLID", current_elem)] +
            //                      (1 - p_phi) * this->elem_mole_fraction_by_phase_[std::make_tuple(
            //                                        id, "LIQUID", current_elem)] -
            //                      current_x;
            // }
          } else {
            mfem::NewtonSolver ourSolver;

            std::unique_ptr<mfem::GMRESSolver> linearSolver = std::make_unique<mfem::GMRESSolver>();
            std::unique_ptr<NonLinearSys_CALPHAD> F = std::make_unique<NonLinearSys_CALPHAD>();

            std::vector<double> lala = {current_x};
            F->get_infos(varphi, lala, c_eq_s, c_eq_l);
            ourSolver.SetOperator(*F);
            ourSolver.SetSolver(*linearSolver);
            ourSolver.SetRelTol(1e-10);
            ourSolver.SetAbsTol(1e-16);
            ourSolver.SetMaxIter(1000);
            ourSolver.SetPrintLevel(0);

            linearSolver->SetRelTol(1e-10);
            linearSolver->SetAbsTol(1e-16);
            linearSolver->SetMaxIter(1000);
            linearSolver->SetPrintLevel(0);

            mfem::Vector x(2);
            x(0) = current_x;
            x(1) = current_x;

            mfem::Vector b(2);
            b = 0.0;

            ourSolver.Mult(b, x);

            // std::cout << "phi : " << varphi << " x_alpha : " << x(0) << " x_beta : " << x(1)
            //           << std::endl;
            // double veri = varphi * x(0) + (1 - varphi) * x(1) - current_x;
            // std::cout << "Verification phi x c_a + (1-phi) x c_b - c =  " << veri << std::endl;

            this->elem_mole_fraction_by_phase_[std::make_tuple(id, "SOLID", current_elem)] = x(0);
            this->elem_mole_fraction_by_phase_[std::make_tuple(id, "LIQUID", current_elem)] = x(1);
          }
          this->chemical_potentials_[std::make_tuple(id, current_elem)] =
              (varphi >= 0.5) ? this->get_mu(this->elem_mole_fraction_by_phase_[std::make_tuple(
                                                id, "SOLID", current_elem)],
                                            c_eq_s, k_solid)
                             : this->get_mu(this->elem_mole_fraction_by_phase_[std::make_tuple(
                                                id, "LIQUID", current_elem)],
                                            c_eq_l, k_liquid);
        }

        // this->elem_mole_fraction_by_phase_[key_mole_fraction_by_phase] =  // current_x;
        // Energy
        for (const auto& energy_name : energy_names) {
          this->energies_of_phases_[std::make_tuple(id, "LIQUID", energy_name)] += this->get_energy(
              this->elem_mole_fraction_by_phase_[std::make_tuple(id, "LIQUID", current_elem)],
              c_eq_l, k_liquid);
          this->energies_of_phases_[std::make_tuple(id, "SOLID", energy_name)] += this->get_energy(
              this->elem_mole_fraction_by_phase_[std::make_tuple(id, "SOLID", current_elem)],
              c_eq_s, k_solid);
        }
      }
    }
    this->driving_force[std::make_tuple(id)] = 0;
    this->driving_force[std::make_tuple(id)] =
        this->energies_of_phases_[std::make_tuple(id, "SOLID", "G")] -
        this->energies_of_phases_[std::make_tuple(id, "LIQUID", "G")];

    for (const auto& element : tuple_x) {
      std::string current_elem = std::get<0>(element);
      this->driving_force[std::make_tuple(id)] -=
          this->chemical_potentials_[std::make_tuple(id, current_elem)] *
          (this->elem_mole_fraction_by_phase_[std::make_tuple(id, "SOLID", current_elem)] -
           this->elem_mole_fraction_by_phase_[std::make_tuple(id, "LIQUID", current_elem)]);
    }
  }
}
/**
 * @brief Check the consistency of outputs required for the current Calphad problem
 *
 * @tparam T
 * @param output_system
 */
template <typename T>
void AnalyticalParaboloidForTwoPhase<T>::check_variables_consistency(
    std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<T>>>& output_system) {
  for (auto& [output_infos, output_value] : output_system) {
    const std::string& output_type = output_infos.back();

    // Fill output with the relevant values
    switch (calphad_outputs::from(output_type)) {
      case calphad_outputs::gm:
      case calphad_outputs::h:
      case calphad_outputs::hm: {
        MFEM_VERIFY(false, "AnalyticalParaboloidForTwoPhase is only built for mu, x and g.\n");

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
void AnalyticalParaboloidForTwoPhase<T>::finalize() {}

/**
 * @brief Destroy the Binary Melting< T>:: Binary Melting object
 *
 * @tparam T
 */
template <typename T>
AnalyticalParaboloidForTwoPhase<T>::~AnalyticalParaboloidForTwoPhase() {}

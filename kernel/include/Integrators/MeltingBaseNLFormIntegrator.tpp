/**
 * @file MeltingBaseNLFormIntegrator.tpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Base class for melting integrators
 * @version 0.1
 * @date 2025-09-05
 *
 * @copyright CEA (C) 2025
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
#include <algorithm>
#include <list>
#include <memory>
#include <span>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Integrators/SlothNLFormIntegrator.hpp"
#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new MeltingBaseNLFormIntegrator object.
 *
 * This constructor initializes the nonlinear form integrator. It forwards the provided previous
 * solution fields, simulation parameters, auxiliary variables, and coefficients to the base SLOTH
 * nonlinear form integrator.
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 *
 * @param u_old        Vector of previous-time-step solution fields.
 * @param params       Paramters that can be used with the integrator.
 * @param auxvars      Auxiliary variables required by the inetgrator.
 * @param coefficients List of coefficients defining material properties.
 *
 */
template <class VARS>
MeltingBaseNLFormIntegrator<VARS>::MeltingBaseNLFormIntegrator(
    const std::vector<mfem::ParGridFunction>& u_old,
    const std::vector<mfem::ParGridFunction>& aux_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients)
    : SlothNLFormIntegrator<VARS>(u_old, aux_old, params, auxvars, coefficients) {}

/**
 * @brief Initialize the integrator.
 *
 * This method performs all necessary setup steps for the integrator:
 * 1. Verifies that the list of expected coefficient types is not empty.
 * 2. Checks that all coefficients contain the expected types.
 * 3. Retrieves and stores the coefficients internally.
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 *
 * @pre The integrator's 'expected_list_' is populated with the
 *      required coefficient types for this integrator.
 *
 */
template <class VARS>
void MeltingBaseNLFormIntegrator<VARS>::init() {
  this->check_coefficient_types(this->expected_list_);
  this->get_coefficients();

  for (std::size_t i = 0; i < this->aux_infos_.size(); ++i) {
    this->vaux_gf_.emplace_back(std::move(this->aux_gf_[i]));
  }
}

/**
 * @brief Assemble the element-level residual vector for the nonlinear problem.
 *
 * This method computes the residual vector  by element
 *
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param el      Array of pointers to finite elements.
 * @param Tr      Element transformation.
 * @param elfun   Array of local finite element solution vectors.
 * @param elvect  Array of vectors where the computed element residual contributions
 *                will be stored.
 *
 * @note Users typically do not call this function directly; it is invoked
 *       internally during the assembly of the global nonlinear form.
 */
template <class VARS>
void MeltingBaseNLFormIntegrator<VARS>::AssembleElementVector(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun, const mfem::Array<mfem::Vector*>& elvect) {
  int num_blocks = el.Size();
  std::vector<double> u_values(2 * num_blocks);

  std::vector<double> vaux_gf_at_ip;
  vaux_gf_at_ip.resize(vaux_gf_.size());
  for (int blk = 0; blk < num_blocks; ++blk) {
    // Catch_Time_Section("MeltingBaseNLFormIntegrator:AssembleElementVector");
    int nd = el[blk]->GetDof();
    int dim = el[blk]->GetDim();
    gradPsi.SetSize(nd, dim);
    Psi.SetSize(nd);
    gradU.SetSize(dim);
    // elvect.SetSize(nd);
    elvect[blk]->SetSize(nd);
    *elvect[blk] = 0.;

    const mfem::IntegrationRule* ir =
        &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());
    // elvect = 0.0;
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      el[blk]->CalcShape(ip, Psi);  //
      Tr.SetIntPoint(&ip);

      const auto& u = *elfun[blk] * Psi;  // Get aux values at ip TODO(cci) (move in method)

      const auto& un = this->u_old_[blk].GetValue(Tr, ip);
      for (size_t k = 0; k < vaux_gf_.size(); ++k) {
        vaux_gf_at_ip[k] = vaux_gf_[k].GetValue(Tr, ip);
      }
      // Get values
      for (int off_blk = 0; off_blk < num_blocks; ++off_blk) {
        u_values[off_blk] = (*elfun[off_blk]) * Psi;
        u_values[off_blk + num_blocks] = this->u_old_[off_blk].GetValue(Tr, ip);
      }

      const double coef_mobi = mobility[blk].compute() * ip.weight * Tr.Weight();
      const double alpha = this->get_phase_change_at_ip(Tr, ip, blk, u, un);
      const double seed = this->get_seed_at_ip(Tr, ip, blk, u, un);
      const double ww = coef_mobi * (alpha * this->compute_gradient_coefficient(
                                                 interpolation_potential[blk], blk,
                                                 std::span<const double>(u_values),
                                                 std::span<const double>(vaux_gf_at_ip)) +
                                     seed);
      add(*elvect[blk], ww, Psi, *elvect[blk]);
    }
  }
}

/**
 * @brief Assemble the element-level Jacobian matrix for the nonlinear problem.
 *
 * This method computes the Jacobian matrix by element.
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param el      Array of pointers to finite elements.
 * @param Tr      Element transformation.
 * @param elfun   Array of local finite element solution vectors.
 * @param elmat   Array of dense matrices where the computed element Jacobian contributions
 *                will be stored.
 *
 * @note Users typically do not call this function directly; it is invoked
 *       internally during the assembly of the global nonlinear form.
 */
template <class VARS>
void MeltingBaseNLFormIntegrator<VARS>::AssembleElementGrad(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun,
    const mfem::Array2D<mfem::DenseMatrix*>& elmats) {
  // Catch_Time_Section("MeltingBaseNLFormIntegrator::AssembleElementGrad");
  int num_blocks = el.Size();
  std::vector<double> u_values(2 * num_blocks);

  std::vector<double> vaux_gf_at_ip;
  vaux_gf_at_ip.resize(vaux_gf_.size());
  for (int blk = 0; blk < num_blocks; ++blk) {
    int nd = el[blk]->GetDof();
    int dim = el[blk]->GetDim();

    gradPsi.SetSize(nd, dim);
    Psi.SetSize(nd);
    // elmat.SetSize(nd);
    elmats(blk, blk)->SetSize(nd);
    *elmats(blk, blk) = 0.0;

    const mfem::IntegrationRule* ir =
        &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());

    // elmat = 0.0;
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      Tr.SetIntPoint(&ip);
      el[blk]->CalcShape(ip, Psi);
      const auto& u = *elfun[blk] * Psi;
      const auto& un = this->u_old_[blk].GetValue(Tr, ip);
      // Get aux values at ip TODO(cci) (move in method)
      for (size_t k = 0; k < vaux_gf_.size(); ++k) {
        vaux_gf_at_ip[k] = vaux_gf_[k].GetValue(Tr, ip);
      }
      // Get values
      for (int off_blk = 0; off_blk < num_blocks; ++off_blk) {
        u_values[off_blk] = (*elfun[off_blk]) * Psi;
        u_values[off_blk + num_blocks] = this->u_old_[off_blk].GetValue(Tr, ip);
      }

      const double coef_mobi = mobility[blk].compute() * ip.weight * Tr.Weight();
      const double alpha = this->get_phase_change_at_ip(Tr, ip, blk, u, un);
      double fun_val =
          coef_mobi * alpha *
          this->compute_hessian_coefficient(
              interpolation_potential[blk], blk, blk, std::span<const double>(u_values),
              std::span<const double>(vaux_gf_at_ip));  // this->energy_derivatives(2, Tr, ip)(u);

      AddMult_a_VVt(fun_val, Psi, *elmats(blk, blk));  // w'(u)*(du, psi)
    }
  }
}

/**
 * @brief Retrieve and store the coefficients required by the integrator.
 *
 * This method collects the coefficients of type  `Mobility`, and `PhaseFieldPotential`
 * from each coefficients and adds them to the corresponding internal storage:
 * - 'mobility' stores the mobility coefficients,
 * - 'interpolation_potential' stores the PhaseFieldPotential coefficients.
 *
 * Only coefficients with ID 0 are considered for each block.
 *
 * @remark This method can be overridden in derived classes to provide
 *         custom behavior for retrieving coefficients.
 *
 * @tparam VARS Template parameter defining the variabls used in the integrator.
 *
 */
template <class VARS>
void MeltingBaseNLFormIntegrator<VARS>::get_coefficients() {
  for (unsigned int i = 0; i < this->nb_blk_; i++) {
    if (this->get_coefficient(i, GlossaryType::Mobility, 0).has_value()) {
      mobility.add(*(this->get_coefficient(i, GlossaryType::Mobility, 0)));
    }
    if (this->get_coefficient(i, GlossaryType::PhaseFieldPotential, 0).has_value()) {
      interpolation_potential.add(
          *(this->get_coefficient(i, GlossaryType::PhaseFieldPotential, 0)));
    }
  }
}

/**
 * @brief Return the value of the coefficient
 * @remark by default values = {u,un} and aux_variables remain accessible in the method with the
 * class variable aux_gf_
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param coef   Coefficient.
 * @param values Vector of current and previous solution values (default: {u, u_old}).

 * @return The computed scalar value of the coefficient.
 */
template <class VARS>
double MeltingBaseNLFormIntegrator<VARS>::compute_coefficient(
    Coefficient coef, const std::span<const double>& values,
    const std::span<const double>& aux_values) {
  if (coef.is_scalar()) {
    return coef.compute();
  } else {
    std::span<const double> u(values.begin(), values.begin() + this->nb_blk_);
    std::span<const double> un(values.begin() + this->nb_blk_, values.end());

    const std::span<const double>& input = coef.is_implicit() ? u : un;
    return coef.compute(input, aux_values);
  }
}

/**
 * @brief Compute the value of a specific component of the gradient of a coefficient.
 *
 * This method evaluates the gradient of the given coefficient with respect to the
 * variable corresponding to a specified index 'blk'. By default, the 'values' vector
 * contains {u, u_old}, and auxiliary variables remain accessible via the class member 'aux_gf_'.
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param coef   Coefficient whose gradient is to be computed.
 * @param blk    Index of the gradient.
 * @param values Vector of current and previous solution values (default: {u, u_old}).
 *
 * @return The computed scalar value of the gradient component of the coefficient.
 *
 * @note This function can be overridden in derived classes.
 */
template <class VARS>
double MeltingBaseNLFormIntegrator<VARS>::compute_gradient_coefficient(
    Coefficient coef, const int blk, const std::span<const double>& values,
    const std::span<const double>& aux_values) {
  std::span<const double> u(values.begin(), values.begin() + this->nb_blk_);
  std::span<const double> un(values.begin() + this->nb_blk_, values.end());
  double coef_value = 0.0;
  if (coef.is_implicit()) {
    coef_value = coef.compute_gradient(blk, u, aux_values);
  } else if (coef.is_explicit()) {
    coef_value = coef.compute_gradient(blk, un, aux_values);
  }
  return coef_value;
}
/**
 * @brief Compute the value of a specific component of the Hessian of a  coefficient.
 *
 * This method evaluates the Hessian of the given coefficient with respect to the
 * variables corresponding to the specified indexes 'iblk' and 'jblk'. By default,
 * the 'values' vector contains {u, u_old}, and auxiliary variables remain accessible
 * via the class member 'aux_gf_'.
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param coef   Coefficient whose Hessian is to be computed.
 * @param iblk   Index of the row  of the Hessian component.
 * @param jblk   Index of the column  of the Hessian component.
 * @param values Vector of current and previous solution values (default: {u, u_old}).
 *
 * @return The computed scalar value of the Hessian component for the given coefficient.
 *
 * @note This function can be overridden in derived classes to implement
 *       more complex or nonlinear Hessian evaluations.
 */
template <class VARS>
double MeltingBaseNLFormIntegrator<VARS>::compute_hessian_coefficient(
    Coefficient coef, const int iblk, const int jblk, const std::span<const double>& values,
    const std::span<const double>& aux_values) {
  std::span<const double> u(values.begin(), values.begin() + this->nb_blk_);
  std::span<const double> un(values.begin() + this->nb_blk_, values.end());
  double coef_value = 0.0;
  if (coef.is_implicit()) {
    coef_value = coef.compute_hessian(iblk, jblk, u, aux_values);
  }
  return coef_value;
}

/**
 * @brief Compute the value of a seed of the secondary phase at integration point
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param Tr Element transformation.
 * @param ir Integration point
 * @param blk   Index of the block.
 * @param u value of the current solution at ip
 * @param un value of the previous solution at ip
 *
 * @return The computed seed at integration point
 */
template <class VARS>
double MeltingBaseNLFormIntegrator<VARS>::get_seed_at_ip(
    [[maybe_unused]] mfem::ElementTransformation& Tr,
    [[maybe_unused]] const mfem::IntegrationPoint& ir, [[maybe_unused]] unsigned int blk,
    [[maybe_unused]] const double u, [[maybe_unused]] const double un) {
  // Nucleus must be equal to zero except when phase transition starts
  const double seed = 0.0;

  return seed;
}

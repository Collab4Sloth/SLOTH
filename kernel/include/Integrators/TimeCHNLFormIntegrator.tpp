/**
 * @file TimeCHNLFormIntegrator.tpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief VF of the time derivative for CH equations
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
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "Coefficients/AxiCylindricalCoefficient.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "Integrators/TimeCHNLFormIntegrator.hpp"
#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new TimeCHNLFormIntegrator object.
 *
 * This constructor initializes a nonlinear form integrator for a splitted time derivative operator.
 * It forwards the provided previous solution fields, simulation
 * parameters, auxiliary variables, and coefficients to the base
 * SLOTH nonlinear form integrator.
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
TimeCHNLFormIntegrator<VARS>::TimeCHNLFormIntegrator(
    Geometry geometry, const double time_step, const std::vector<mfem::ParGridFunction>& u_old,
    const std::vector<mfem::ParGridFunction>& aux_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients)
    : SlothNLFormIntegrator<VARS>(geometry, time_step, u_old, aux_old, params, auxvars,
                                  coefficients) {
  this->integrator_name_ = "SplitTimeDerivative";

  this->check_variables_consistency();
}

/**
 * @brief Initialize the SplitTimeDerivative integrator.
 *
 * This method performs all necessary setup steps for the SplitTimeDerivative
 * nonlinear form integrator:
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
void TimeCHNLFormIntegrator<VARS>::init() {
  if (this->expected_list_.empty()) {
    this->get_coefficients();
  } else {
    this->check_coefficient_types(this->expected_list_);
    this->get_coefficients();
  }

  for (std::size_t i = 0; i < this->aux_infos_.size(); ++i) {
    this->vaux_gf_.emplace_back(std::move(this->aux_gf_[i]));
  }
}

/**
 * @brief  Check variables consistency
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 */
template <class VARS>
void TimeCHNLFormIntegrator<VARS>::check_variables_consistency() {
  // Temperature scaling for mobility
  for (std::size_t i = 0; i < this->aux_infos_.size(); ++i) {
    const auto& variable_info = this->aux_infos_[i];
    MFEM_VERIFY(!variable_info.empty(), "Empty variable_info encountered.");
    size_t vsize = variable_info.size();

    MFEM_VERIFY(vsize >= 1,
                "TimeCHNLFormIntegrator<VARS>: at least "
                "one additionnal information is expected for auxiliary variables associated with "
                "this integrator");
    const std::string& symbol = toUpperCase(variable_info.back());
    if (symbol == "T") {
      this->temp_gf_.emplace_back(std::move(this->aux_gf_[i]));
      break;
    }
  }
}

/**
 * @brief Retrieve and store the coefficients required by the TimeDerivative integrator.
 *
 * This method collects one default coefficient:
 * - 'coefficient_A' stores the first coefficient,
 *
 *
 * @remark By default, A is equal to one. This method can be overridden in derived classes to
 * provide custom behavior for retrieving coefficients. Coefficient object
 * of type GlossaryType::ExplicitTime can be also used for this class (see Glossary::ExplicitTime_A
 * used wth explicit solver)
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 */
template <class VARS>
void TimeCHNLFormIntegrator<VARS>::get_coefficients() {
  for (unsigned int i = 0; i < this->nb_blk_; i++) {
    auto coef_A = this->get_coefficient(i, GlossaryType::ExplicitTime, 0);
    this->coefficient_A.add(coef_A.has_value() ? *coef_A : Coefficient(Glossary::Default, 1.0));
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
void TimeCHNLFormIntegrator<VARS>::AssembleElementVector(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun, const mfem::Array<mfem::Vector*>& elvect) {
  // {u, u_old} and auxiliary fields at the integration points (same rule as below)
  this->prepare_element_values(
      el, Tr, elfun, this->vaux_gf_,
      mfem::IntRules.Get(el[0]->GetGeomType(), 2 * el[0]->GetOrder() + Tr.OrderW()));
  //////////////////////
  // Block 0 R(phi) on mu term
  {
    int blk = 0;

    int nd = el[blk]->GetDof();
    elvect[blk]->SetSize(nd);
    *elvect[blk] = 0.;
  }

  //////////////////////
  // Block 1 R(mu) on dphi/dt term
  {
    int blk = 1;

    int off_blk = 0;

    int nd = el[blk]->GetDof();
    int dim = el[blk]->GetDim();
    gradPsi.SetSize(nd, dim);
    Psi.SetSize(nd);
    gradU.SetSize(dim);
    elvect[blk]->SetSize(nd);
    *elvect[blk] = 0.;

    const mfem::IntegrationRule* ir =
        &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      el[blk]->CalcShape(ip, Psi);  //
      this->set_int_point(Tr, ip);

      const auto& phi = *elfun[off_blk] * Psi;
      const auto ipv = this->values_at_ip(*ir, i, Psi, elfun, this->vaux_gf_, Tr);

      double weight_coef = ip.weight * Tr.Weight();
      if (this->isAxisymmetric()) {
        weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
      }
      // Implicit coefficient: evaluated at the end of the step (phi = phi_old + dt dphi/dt)
      const double coef_a =
          this->eval_coefficient(this->coefficient_A[blk], this->values_at_end_of_step(ipv.u),
                                 ipv.aux);
      const double ww = coef_a * phi * weight_coef;
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
void TimeCHNLFormIntegrator<VARS>::AssembleElementGrad(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun,
    const mfem::Array2D<mfem::DenseMatrix*>& elmats) {
  // {k, u_old} and auxiliary fields at the integration points (same rule as below)
  this->prepare_element_values(
      el, Tr, elfun, this->vaux_gf_,
      mfem::IntRules.Get(el[0]->GetGeomType(), 2 * el[0]->GetOrder() + Tr.OrderW()));
  for (int blk = 0; blk < 2; ++blk) {
    for (int jblk = 0; jblk < 2; ++jblk) {
      elmats(blk, jblk)->SetSize(el[blk]->GetDof(), el[jblk]->GetDof());
      *elmats(blk, jblk) = 0.0;
    }
  }

  // Row 1 (equation of mu), term (A(u) dphi/dt, psi):
  //   block (1, 0): A (psi_j, psi_i) + dt dA/dphi dphi/dt (psi_j, psi_i)
  //   block (1, 1): dt dA/dmu dphi/dt (psi_j, psi_i)
  // (u = u_old + dt k at the end of the step; the derivatives vanish for constant or explicit A)
  const int blk = 1;
  const int off_blk = 0;
  Psi.SetSize(el[blk]->GetDof());
  const mfem::IntegrationRule* ir =
      &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());
  for (int i = 0; i < ir->GetNPoints(); i++) {
    const mfem::IntegrationPoint& ip = ir->IntPoint(i);
    this->set_int_point(Tr, ip);
    el[blk]->CalcShape(ip, Psi);

    const auto ipv = this->values_at_ip(*ir, i, Psi, elfun, this->vaux_gf_, Tr);
    const auto u_end = this->values_at_end_of_step(ipv.u);

    double weight_coef = ip.weight * Tr.Weight();
    if (this->isAxisymmetric()) {
      weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
    }

    const double coef_a = this->eval_coefficient(this->coefficient_A[blk], u_end, ipv.aux);
    AddMult_a_VVt(coef_a * weight_coef, Psi, *elmats(blk, off_blk));

    const double dphi_dt = *elfun[off_blk] * Psi;
    for (int jblk = 0; jblk < 2; ++jblk) {
      const double d =
          this->eval_gradient_coefficient(this->coefficient_A[blk], jblk, u_end, ipv.aux) *
          this->time_step_ * dphi_dt * weight_coef;
      if (d != 0.0) AddMult_a_VVt(d, Psi, *elmats(blk, jblk));
    }
  }
}

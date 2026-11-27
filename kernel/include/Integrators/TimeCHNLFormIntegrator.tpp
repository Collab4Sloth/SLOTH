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

  // A of the equations of mu_i (blocks n + i, coefficient of one variable: c_i) evaluated at the
  // end of the step at all the quadrature points of the local mesh before each assembly, with the
  // outputs of the residual and of the Jacobian
  const unsigned int value = SlothNLBase::CoefValue;
  const unsigned int first = SlothNLBase::CoefFirstDerivative;
  const int n = static_cast<int>(this->nb_blk_) / 2;
  for (int off_blk = 0; off_blk < n; ++off_blk) {
    this->prepare_on_mesh(this->coefficient_A[n + off_blk], off_blk, value, value | first,
                          this->vaux_gf_, SlothNLBase::ValuesAt::EndOfStep);
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
  // {u, u_old} and auxiliary fields at the integration points (same rule as below): not
  // needed if the coefficients are prepared on the mesh
  const mfem::IntegrationRule& ir_values = this->integration_rule(*el[0], Tr);
  if (!this->values_prepared_on_mesh(Tr.ElementNo, ir_values)) {
    this->prepare_element_values(el, Tr, elfun, this->vaux_gf_, ir_values);
  }
  // Blocks c_0..c_{n-1}, mu_0..mu_{n-1} (n = 1: phi, mu)
  const int n = el.Size() / 2;
  //////////////////////
  // Blocks c_i: R(c_i) on mu_i term
  for (int blk = 0; blk < n; ++blk) {
    int nd = el[blk]->GetDof();
    elvect[blk]->SetSize(nd);
    *elvect[blk] = 0.;
  }

  //////////////////////
  // Blocks mu_i: R(mu_i) on dc_i/dt term
  for (int off_blk = 0; off_blk < n; ++off_blk) {
    int blk = n + off_blk;

    int nd = el[blk]->GetDof();
    int dim = el[blk]->GetDim();
    gradPsi.SetSize(nd, dim);
    Psi.SetSize(nd);
    gradU.SetSize(dim);
    elvect[blk]->SetSize(nd);
    *elvect[blk] = 0.;

    const mfem::IntegrationRule* ir = &this->integration_rule(*el[blk], Tr);
    // Implicit coefficient: evaluated at the end of the step (c = c_old + dt dc/dt), at all the
    // integration points of the element; a coefficient of one variable is a function of c_i
    this->eval_coefficient_at_points(this->coefficient_A[blk], off_blk, SlothNLBase::CoefValue,
                                     *el[blk], Tr, *ir, elfun, this->vaux_gf_, this->a_at_points_,
                                     SlothNLBase::ValuesAt::EndOfStep);
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      el[blk]->CalcShape(ip, Psi);  //
      this->set_int_point(Tr, ip);

      const auto& phi = *elfun[off_blk] * Psi;

      double weight_coef = ip.weight * Tr.Weight();
      if (this->isAxisymmetric()) {
        weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
      }
      const double coef_a = this->a_at_points_.value(i);
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
  // {k, u_old} and auxiliary fields at the integration points (same rule as below): not
  // needed if the coefficients are prepared on the mesh
  const mfem::IntegrationRule& ir_values = this->integration_rule(*el[0], Tr);
  if (!this->values_prepared_on_mesh(Tr.ElementNo, ir_values)) {
    this->prepare_element_values(el, Tr, elfun, this->vaux_gf_, ir_values);
  }
  // Blocks c_0..c_{n-1}, mu_0..mu_{n-1} (n = 1: phi, mu)
  const int nb = el.Size();
  const int n = nb / 2;
  for (int blk = 0; blk < nb; ++blk) {
    for (int jblk = 0; jblk < nb; ++jblk) {
      elmats(blk, jblk)->SetSize(el[blk]->GetDof(), el[jblk]->GetDof());
      *elmats(blk, jblk) = 0.0;
    }
  }

  // Row n + i (equation of mu_i), term (A(u) dc_i/dt, psi):
  //   block (n + i, i): A (psi_j, psi_i) + dt dA/dc_i dc_i/dt (psi_j, psi_i)
  //   block (n + i, j): dt dA/du_j dc_i/dt (psi_j, psi_i)
  // (u = u_old + dt k at the end of the step; the derivatives vanish for constant or explicit A)
  for (int off_blk = 0; off_blk < n; ++off_blk) {
    const int blk = n + off_blk;
    Psi.SetSize(el[blk]->GetDof());
    const mfem::IntegrationRule* ir = &this->integration_rule(*el[blk], Tr);
    // A and its derivatives at all the integration points of the element, at the end of the step
    this->eval_coefficient_at_points(this->coefficient_A[blk], off_blk,
                                     SlothNLBase::CoefValue | SlothNLBase::CoefFirstDerivative,
                                     *el[blk], Tr, *ir, elfun, this->vaux_gf_, this->a_at_points_,
                                     SlothNLBase::ValuesAt::EndOfStep);
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      this->set_int_point(Tr, ip);
      el[blk]->CalcShape(ip, Psi);

      double weight_coef = ip.weight * Tr.Weight();
      if (this->isAxisymmetric()) {
        weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
      }

      const double coef_a = this->a_at_points_.value(i);
      AddMult_a_VVt(coef_a * weight_coef, Psi, *elmats(blk, off_blk));

      const double dphi_dt = *elfun[off_blk] * Psi;
      for (int jblk = 0; jblk < nb; ++jblk) {
        const double d =
            this->a_at_points_.derivative(i, jblk) * this->time_step_ * dphi_dt * weight_coef;
        if (d != 0.0) AddMult_a_VVt(d, Psi, *elmats(blk, jblk));
      }
    }
  }
}

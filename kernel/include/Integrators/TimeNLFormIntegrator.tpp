/**
 * @file TimeNLFormIntegrator.tpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief VF of the time derivative
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
#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new TimeNLFormIntegrator object.
 *
 * This constructor initializes a nonlinear form integrator for a simple time derivative operator.
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
TimeNLFormIntegrator<VARS>::TimeNLFormIntegrator(Geometry geometry, const double time_step,
                                                 const std::vector<mfem::ParGridFunction> u_old,
                                                 const std::vector<mfem::ParGridFunction> aux_old,
                                                 const Parameters& params,
                                                 std::vector<VARS*> auxvars,
                                                 const std::vector<Coefficients>& coefficients)
    : SlothNLFormIntegrator<VARS>(geometry, time_step, u_old, aux_old, params, auxvars,
                                  coefficients) {
  this->integrator_name_ = "TimeDerivative";

  this->check_variables_consistency();
}

/**
 * @brief Initialize the TimeDerivative integrator.
 *
 * This method performs all necessary setup steps for the TimeDerivative
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
void TimeNLFormIntegrator<VARS>::init() {
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
void TimeNLFormIntegrator<VARS>::check_variables_consistency() {
  // Temperature scaling for mobility
  for (std::size_t i = 0; i < this->aux_infos_.size(); ++i) {
    const auto& variable_info = this->aux_infos_[i];
    MFEM_VERIFY(!variable_info.empty(), "Empty variable_info encountered.");
    size_t vsize = variable_info.size();

    MFEM_VERIFY(vsize >= 1,
                "TimeNLFormIntegrator<VARS>: at least "
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
 * This method collects two default coefficients:
 * - 'coefficient_A' stores the first coefficient,
 * - 'coefficient_B' stores the second coefficient,
 *
 *
 * @remark By default, A and B are equal to one. Implicit coefficients are evaluated at the end
 * of the time step (see SlothNLFormIntegrator::values_at_end_of_step). This method can be
 * overridden in derived classes to
 * provide custom behavior for retrieving coefficients (see HeatTimeDerivative). Coefficient objects
 * of type GlossaryType::ExplicitTime can be also used for this class (see Glossary::ExplicitTime_A
 * and Glossary::ExplicitTime_B used wth explicit solver)
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 */
template <class VARS>
void TimeNLFormIntegrator<VARS>::get_coefficients() {
  for (unsigned int i = 0; i < this->nb_blk_; i++) {
    auto coef_A = this->get_coefficient(i, GlossaryType::ExplicitTime, 0);
    this->coefficient_A.add(coef_A.has_value() ? *coef_A : Coefficient(Glossary::Default, 1.0));
    auto coef_B = this->get_coefficient(i, GlossaryType::ExplicitTime, 1);
    this->coefficient_B.add(coef_B.has_value() ? *coef_B : Coefficient(Glossary::Default, 1.0));
  }
}

/**
 * @brief Element residual: (A(u) B(u) k, psi) for each block, with k = du/dt the unknown and
 * u = u_old + dt k the value at the end of the step (A and B: implicit, explicit or constant
 * coefficients).
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 */
template <class VARS>
void TimeNLFormIntegrator<VARS>::AssembleElementVector(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun, const mfem::Array<mfem::Vector*>& elvect) {
  const int num_blocks = el.Size();
  // {k, u_old} and auxiliary fields at the integration points (same rule as below)
  this->prepare_element_values(
      el, Tr, elfun, this->vaux_gf_,
      mfem::IntRules.Get(el[0]->GetGeomType(), 2 * el[0]->GetOrder() + Tr.OrderW()));
  for (int blk = 0; blk < num_blocks; ++blk) {
    const int nd = el[blk]->GetDof();
    Psi.SetSize(nd);
    elvect[blk]->SetSize(nd);
    *elvect[blk] = 0.;

    const mfem::IntegrationRule* ir =
        &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      el[blk]->CalcShape(ip, Psi);
      this->set_int_point(Tr, ip);

      const double k = *elfun[blk] * Psi;
      const auto ipv = this->values_at_ip(*ir, i, Psi, elfun, this->vaux_gf_, Tr);
      const auto u_end = this->values_at_end_of_step(ipv.u);
      const double coef_a = this->eval_coefficient(this->coefficient_A[blk], u_end, ipv.aux, blk);
      const double coef_b = this->eval_coefficient(this->coefficient_B[blk], u_end, ipv.aux, blk);

      double weight_coef = ip.weight * Tr.Weight();
      if (this->isAxisymmetric()) {
        weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
      }
      elvect[blk]->Add(coef_a * coef_b * k * weight_coef, Psi);
    }
  }
}

/**
 * @brief Element Jacobian with respect to k:
 *   block (b, b): A B (psi_j, psi_i) + dt d(A B)/du_b k (psi_j, psi_i)
 *   block (b, j): dt d(A B)/du_j k (psi_j, psi_i)
 * since du/dk = dt (u = u_old + dt k). The derivatives vanish for constant or explicit
 * coefficients.
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 */
template <class VARS>
void TimeNLFormIntegrator<VARS>::AssembleElementGrad(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun,
    const mfem::Array2D<mfem::DenseMatrix*>& elmats) {
  const int num_blocks = el.Size();
  // {k, u_old} and auxiliary fields at the integration points (same rule as below)
  this->prepare_element_values(
      el, Tr, elfun, this->vaux_gf_,
      mfem::IntRules.Get(el[0]->GetGeomType(), 2 * el[0]->GetOrder() + Tr.OrderW()));
  for (int blk = 0; blk < num_blocks; ++blk) {
    for (int jblk = 0; jblk < num_blocks; ++jblk) {
      elmats(blk, jblk)->SetSize(el[blk]->GetDof(), el[jblk]->GetDof());
      *elmats(blk, jblk) = 0.0;
    }
  }

  mfem::Vector shape_j;
  for (int blk = 0; blk < num_blocks; ++blk) {
    Psi.SetSize(el[blk]->GetDof());
    const mfem::IntegrationRule* ir =
        &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      el[blk]->CalcShape(ip, Psi);
      this->set_int_point(Tr, ip);

      const double k = *elfun[blk] * Psi;
      const auto ipv = this->values_at_ip(*ir, i, Psi, elfun, this->vaux_gf_, Tr);
      const auto u_end = this->values_at_end_of_step(ipv.u);
      const double coef_a = this->eval_coefficient(this->coefficient_A[blk], u_end, ipv.aux, blk);
      const double coef_b = this->eval_coefficient(this->coefficient_B[blk], u_end, ipv.aux, blk);

      double weight_coef = ip.weight * Tr.Weight();
      if (this->isAxisymmetric()) {
        weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
      }
      AddMult_a_VVt(coef_a * coef_b * weight_coef, Psi, *elmats(blk, blk));

      // Coefficients depending on u at the end of the step: d/dk_j = dt d/du_j
      for (int jblk = 0; jblk < num_blocks; ++jblk) {
        const double d_a =
            this->eval_gradient_coefficient(this->coefficient_A[blk], jblk, u_end, ipv.aux, blk);
        const double d_b =
            this->eval_gradient_coefficient(this->coefficient_B[blk], jblk, u_end, ipv.aux, blk);
        const double d = (d_a * coef_b + coef_a * d_b) * this->time_step_ * k * weight_coef;
        if (d == 0.0) continue;
        if (jblk == blk) {
          AddMult_a_VVt(d, Psi, *elmats(blk, blk));
        } else {
          shape_j.SetSize(el[jblk]->GetDof());
          el[jblk]->CalcShape(ip, shape_j);
          AddMult_a_VWt(d, Psi, shape_j, *elmats(blk, jblk));
        }
      }
    }
  }
}

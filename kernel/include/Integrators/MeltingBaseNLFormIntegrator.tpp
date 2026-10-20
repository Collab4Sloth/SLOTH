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

#include "Coefficients/AxiCylindricalCoefficient.hpp"
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
    Geometry geometry, const double time_step, const std::vector<mfem::ParGridFunction>& u_old,
    const std::vector<mfem::ParGridFunction>& aux_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients)
    : SlothNLFormIntegrator<VARS>(geometry, time_step, u_old, aux_old, params, auxvars,
                                  coefficients) {}

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
 * @brief  Check variables consistency
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 */
template <class VARS>
void MeltingBaseNLFormIntegrator<VARS>::check_variables_consistency() {}

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
  // {u, u_old} and auxiliary fields at the integration points (same rule as below)
  this->prepare_element_values(
      el, Tr, elfun, this->vaux_gf_,
      mfem::IntRules.Get(el[0]->GetGeomType(), 2 * el[0]->GetOrder() + Tr.OrderW()));
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
      this->set_int_point(Tr, ip);

      const auto ipv = this->values_at_ip(*ir, i, Psi, elfun, this->vaux_gf_, Tr);

      const double coef_mobi = get_mob_at_ip(blk, ipv.u, ipv.aux);
      const double alpha = this->get_phase_change_at_ip(blk, ipv.u, ipv.aux);
      const double seed = this->get_seed_at_ip(blk, ipv.u, ipv.aux);
      double ww = coef_mobi * (alpha * this->eval_gradient_energy_coefficient(
                                           interpolation_potential[blk], blk, ipv.u, ipv.aux)) +
                  seed;

      double weight_coef = ip.weight * Tr.Weight();
      if (this->isAxisymmetric()) {
        weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
      }

      ww *= weight_coef;
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
  // {u, u_old} and auxiliary fields at the integration points (same rule as below)
  this->prepare_element_values(
      el, Tr, elfun, this->vaux_gf_,
      mfem::IntRules.Get(el[0]->GetGeomType(), 2 * el[0]->GetOrder() + Tr.OrderW()));
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
      this->set_int_point(Tr, ip);
      el[blk]->CalcShape(ip, Psi);
      const auto ipv = this->values_at_ip(*ir, i, Psi, elfun, this->vaux_gf_, Tr);

      const double coef_mobi = get_mob_at_ip(blk, ipv.u, ipv.aux);
      const double alpha = this->get_phase_change_at_ip(blk, ipv.u, ipv.aux);
      double fun_val =
          coef_mobi * alpha *
          this->eval_hessian_coefficient(interpolation_potential[blk], blk, blk, ipv.u,
                                         ipv.aux);  // this->energy_derivatives(2, Tr, ip)(u);

      double weight_coef = ip.weight * Tr.Weight();
      if (this->isAxisymmetric()) {
        weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
      }

      fun_val *= weight_coef;

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
    [[maybe_unused]] unsigned int blk, [[maybe_unused]] const std::span<const double>& values,
    [[maybe_unused]] const std::span<const double>& aux_values) {
  // Nucleus must be equal to zero except when phase transition starts
  const double seed = 0.0;

  return seed;
}

/**
 * @brief Compute the value of the mobility at integration point
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param Tr Element transformation.
 * @param ir Integration point
 * @param blk   Index of the block.
 * @param u value of the current solution at ip
 * @param un value of the previous solution at ip
 *
 * @return The computed mobility at integration point
 */
template <class VARS>
double MeltingBaseNLFormIntegrator<VARS>::get_mob_at_ip(
    [[maybe_unused]] unsigned int blk, [[maybe_unused]] const std::span<const double>& values,
    [[maybe_unused]] const std::span<const double>& aux_values) {
  const double mobi = this->eval_coefficient(mobility[blk], values, aux_values);

  return mobi;
}

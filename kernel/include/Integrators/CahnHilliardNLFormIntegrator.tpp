/**
 * @file CahnHilliardNLFormIntegrator.tpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief VF of the CahnHilliard equations
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
#include <string>
#include <utility>
#include <vector>

#include "Coefficients/AxiCylindricalCoefficient.hpp"
#include "Coefficients/SlothBaseCoefficient.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new CahnHilliardNLFormIntegrator object.
 *
 * This constructor initializes a nonlinear form integrator for Cahn-Hilliard-type
 * problems. It forwards the provided previous solution fields, simulation
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
CahnHilliardNLFormIntegrator<VARS>::CahnHilliardNLFormIntegrator(
    Geometry geometry, const double time_step, const std::vector<mfem::ParGridFunction>& u_old,
    const std::vector<mfem::ParGridFunction>& aux_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients)
    : SlothNLFormIntegrator<VARS>(geometry, time_step, u_old, aux_old, params, auxvars,
                                  coefficients) {
  this->integrator_name_ = "CahnHilliard";
  this->check_variables_consistency();
}

/**
 * @brief Initialize the CahnHilliard integrator.
 *
 * This method performs all necessary setup steps for a CahnHilliard-type
 * nonlinear form integrator:
 * 1. Verifies that the list of expected coefficient types is not empty.
 * 2. Checks that all coefficients contain the expected types.
 * 3. Retrieves and stores the coefficients internally.
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 *
 * @pre The integrator's 'expected_list_' is populated with the
 *      required coefficient types for this CahnHilliard integrator.
 *
 */
template <class VARS>
void CahnHilliardNLFormIntegrator<VARS>::init() {
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
void CahnHilliardNLFormIntegrator<VARS>::check_variables_consistency() {
  // TODO(cci): order of variables?
}

/**
 * @brief Retrieve and store the coefficients required by the Cahn-Hilliard integrator.
 *
 * This method collects the coefficients of type `Capillary`, `Mobility`, and `FreeEnergy`
 * from each coefficients and adds them to the corresponding internal storage:
 * - 'lambda' stores the capillary coefficients,
 * - 'mobility' stores the mobility coefficients,
 * - 'double_well_energy' stores the free energy coefficients.
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
void CahnHilliardNLFormIntegrator<VARS>::get_coefficients() {
  for (unsigned int i = 0; i < this->nb_blk_; i++) {
    if (this->get_coefficient(i, GlossaryType::Capillary, 0).has_value()) {
      lambda.add(*(this->get_coefficient(i, GlossaryType::Capillary, 0)));
    }
    if (this->get_coefficient(i, GlossaryType::Mobility, 0).has_value()) {
      mobility.add(*(this->get_coefficient(i, GlossaryType::Mobility, 0)));
    }
    if (this->get_coefficient(i, GlossaryType::FreeEnergy, 0).has_value()) {
      double_well_energy.add(*(this->get_coefficient(i, GlossaryType::FreeEnergy, 0)));
    }
  }
}

/**
 * @brief Assemble the element-level residual vector for the nonlinear problem.
 *
 * This method computes the residual vector of the CahnHilliard-type nonlinear problem by element
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
void CahnHilliardNLFormIntegrator<VARS>::AssembleElementVector(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun, const mfem::Array<mfem::Vector*>& elvect) {
  // {u, u_old} and auxiliary fields at the integration points (same rule as below)
  this->prepare_element_values(
      el, Tr, elfun, this->vaux_gf_,
      mfem::IntRules.Get(el[0]->GetGeomType(), 2 * el[0]->GetOrder() + Tr.OrderW()));

  //
  // Block 0 R(phi)=mu - w' + div lambda grad phi
  //
  {
    int blk = 0;
    int off_blk = 1;
    mfem::DenseMatrix gradPsi;
    mfem::Vector Psi, gradU;
    int nd = el[blk]->GetDof();
    int dim = el[blk]->GetDim();
    gradPsi.SetSize(nd, dim);
    Psi.SetSize(nd);
    gradU.SetSize(dim);
    elvect[blk]->SetSize(nd);
    *elvect[blk] = 0.;
    const mfem::IntegrationRule* ir =
        &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());
    ///// CCI
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      el[blk]->CalcShape(ip, Psi);  //
      this->set_int_point(Tr, ip);

      const auto& mu = *elfun[off_blk] * Psi;
      const auto ipv = this->values_at_ip(*ir, i, Psi, elfun, this->vaux_gf_, Tr);

      double weight_coef = ip.weight * Tr.Weight();
      if (this->isAxisymmetric()) {
        weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
      }
      el[blk]->CalcPhysDShape(Tr, gradPsi);
      gradPsi.MultTranspose(*elfun[blk], gradU);

      gradU *= -weight_coef * this->eval_coefficient(lambda[blk], ipv.u, ipv.aux);
      gradPsi.AddMult(gradU, *elvect[blk]);

      // Given u, compute (w'(u), psi), psi is shape function
      const double ww = weight_coef * (mu - this->eval_gradient_energy_coefficient(
                                                double_well_energy[blk], blk, ipv.u, ipv.aux));

      add(*elvect[blk], ww, Psi, *elvect[blk]);
    }
  }
  //
  // Block 1 R(mu) = div M grad mu
  //
  {
    int blk = 1;
    mfem::DenseMatrix gradPsi;
    mfem::Vector Psi, gradU;
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

      const auto ipv = this->values_at_ip(*ir, i, Psi, elfun, this->vaux_gf_, Tr);

      double weight_coef = ip.weight * Tr.Weight();
      if (this->isAxisymmetric()) {
        weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
      }

      el[blk]->CalcPhysDShape(Tr, gradPsi);
      gradPsi.MultTranspose(*elfun[blk], gradU);
      gradU *= this->eval_coefficient(mobility[blk], ipv.u, ipv.aux) * weight_coef;
      gradPsi.AddMult(gradU, *elvect[blk]);
    }
  }
}

/**
 * @brief Assemble the element-level Jacobian matrix for the nonlinear problem.
 *
 * This method computes the Jacobian matrix of a CahHilliard-type nonlinear problem by element.
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
void CahnHilliardNLFormIntegrator<VARS>::AssembleElementGrad(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun,
    const mfem::Array2D<mfem::DenseMatrix*>& elmats) {
  // {u, u_old} and auxiliary fields at the integration points (same rule as below)
  this->prepare_element_values(
      el, Tr, elfun, this->vaux_gf_,
      mfem::IntRules.Get(el[0]->GetGeomType(), 2 * el[0]->GetOrder() + Tr.OrderW()));

  // Block 0  0 dR(phi)dphi = d(mu - w' + div lambda grad phi)/dphi
  {
    int blk = 0;
    mfem::DenseMatrix gradPsi;
    mfem::Vector Psi, gradU;
    int nd = el[blk]->GetDof();
    int dim = el[blk]->GetDim();

    gradPsi.SetSize(nd, dim);
    Psi.SetSize(nd);

    elmats(blk, blk)->SetSize(nd);
    *elmats(blk, blk) = 0.0;

    const mfem::IntegrationRule* ir =
        &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      el[blk]->CalcShape(ip, Psi);
      this->set_int_point(Tr, ip);

      const auto ipv = this->values_at_ip(*ir, i, Psi, elfun, this->vaux_gf_, Tr);

      double weight_coef = ip.weight * Tr.Weight();
      if (this->isAxisymmetric()) {
        weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
      }

      el[blk]->CalcPhysDShape(Tr, gradPsi);

      const double coef_lambda = this->eval_coefficient(lambda[blk], ipv.u, ipv.aux);

      AddMult_a_AAt(-weight_coef * coef_lambda, gradPsi, *elmats(blk, blk));

      double fun_val = -weight_coef * this->eval_hessian_coefficient(double_well_energy[blk], blk,
                                                                     blk, ipv.u, ipv.aux);
      AddMult_a_VVt(fun_val, Psi, *elmats(blk, blk));
    }
  }

  // Block 0 1 dR(phi)dmu=d(mu - w' + div lambda grad phi)/dmu
  {
    int blk = 0;
    int off_blk = 1;
    mfem::DenseMatrix gradPsi;
    mfem::Vector Psi, gradU;
    int nd = el[blk]->GetDof();
    int dim = el[blk]->GetDim();

    gradPsi.SetSize(nd, dim);
    Psi.SetSize(nd);

    elmats(blk, off_blk)->SetSize(nd);
    *elmats(blk, off_blk) = 0.0;
    const mfem::IntegrationRule* ir =
        &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      el[blk]->CalcShape(ip, Psi);
      this->set_int_point(Tr, ip);

      double weight_coef = ip.weight * Tr.Weight();
      if (this->isAxisymmetric()) {
        weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
      }

      AddMult_a_VVt(weight_coef, Psi, *elmats(blk, off_blk));
    }
  }
  // Block 1 0  dR(mu)dphi=d(div M grad mu)/dphi
  {
    int blk = 1;
    int off_blk = 0;
    mfem::DenseMatrix gradPsi;
    mfem::Vector Psi, gradU;
    int nd = el[blk]->GetDof();
    int dim = el[blk]->GetDim();

    gradPsi.SetSize(nd, dim);
    Psi.SetSize(nd);

    elmats(blk, off_blk)->SetSize(nd, el[off_blk]->GetDof());
    *elmats(blk, off_blk) = 0.0;

    // Mobility depending on phi:
    //   dR(mu)_i/dphi_j = (dM/dphi) psi_j grad(mu) . grad(psi_i)
    // (zero, and skipped, for a mobility independent of phi)
    mfem::Vector Psi_phi(el[off_blk]->GetDof());
    mfem::Vector grad_mu(dim), adv(nd);
    const mfem::IntegrationRule* ir =
        &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      el[blk]->CalcShape(ip, Psi);
      this->set_int_point(Tr, ip);
      const auto ipv = this->values_at_ip(*ir, i, Psi, elfun, this->vaux_gf_, Tr);

      double weight_coef = ip.weight * Tr.Weight();
      if (this->isAxisymmetric()) {
        weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
      }
      const double dmob_dphi =
          this->eval_gradient_coefficient(mobility[blk], off_blk, ipv.u, ipv.aux) * weight_coef;
      if (dmob_dphi == 0.0) continue;

      el[blk]->CalcPhysDShape(Tr, gradPsi);
      gradPsi.MultTranspose(*elfun[blk], grad_mu);
      gradPsi.Mult(grad_mu, adv);  // adv_i = grad(mu) . grad(psi_i)
      el[off_blk]->CalcShape(ip, Psi_phi);
      AddMult_a_VWt(dmob_dphi, adv, Psi_phi, *elmats(blk, off_blk));
    }
  }

  // Block 1 1  dR(mu)dmu=dR(mu)dphi=d(div M grad mu)/dmu
  {
    int blk = 1;
    mfem::DenseMatrix gradPsi;
    mfem::Vector Psi, gradU;
    int nd = el[blk]->GetDof();
    int dim = el[blk]->GetDim();

    gradPsi.SetSize(nd, dim);
    Psi.SetSize(nd);

    elmats(blk, blk)->SetSize(nd);
    *elmats(blk, blk) = 0.0;
    const mfem::IntegrationRule* ir =
        &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());
    ///// CCI
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      el[blk]->CalcShape(ip, Psi);  //
      this->set_int_point(Tr, ip);

      const auto ipv = this->values_at_ip(*ir, i, Psi, elfun, this->vaux_gf_, Tr);

      double weight_coef = ip.weight * Tr.Weight();
      if (this->isAxisymmetric()) {
        weight_coef *= AxiCylindricalCoefficient().Eval(Tr, ip);
      }

      const double coef_mob = this->eval_coefficient(mobility[blk], ipv.u, ipv.aux) * weight_coef;
      el[blk]->CalcPhysDShape(Tr, gradPsi);

      AddMult_a_AAt(coef_mob, gradPsi, *elmats(blk, blk));
    }
  }
}

/**
 * @file DiffusionFluxNLFormIntegrator.tpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief inter-diffusion integrator
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

#include <list>
#include <memory>
#include <span>
#include <utility>
#include <vector>

#include "Coefficients/SlothBaseCoefficient.hpp"
#include "Integrators/SlothGridFunction.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameters.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Retrieve and store the coefficients required by the Allen-Cahn integrator.
 *
 * This method collects a coefficient of type Diffusivity and stores it:
 * - 'stab_diffusion' stores the coefficient use to stabilize the inter-diffusion equation.
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
void DiffusionFluxNLFormIntegrator<VARS>::get_coefficients() {
  for (unsigned int i = 0; i < this->nb_blk_; i++) {
    if (this->get_coefficient(i, GlossaryType::Diffusivity, 0).has_value()) {
      this->stab_diffusion.add(*(this->get_coefficient(i, GlossaryType::Diffusivity, 0)));
    }
  }
}

/**
 * @brief Construct a new DiffusionFluxNLFormIntegrator object.
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
DiffusionFluxNLFormIntegrator<VARS>::DiffusionFluxNLFormIntegrator(
    Geometry geometry, const std::vector<mfem::ParGridFunction>& u_old,
    const std::vector<mfem::ParGridFunction>& aux_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients)
    : SlothNLFormIntegrator<VARS>(geometry, u_old, aux_old, params, auxvars, coefficients) {
  this->expected_list_.push_back(GlossaryType::Diffusivity);
}

/**
 * @brief Initialize the integrator.
 *
 * This method performs all necessary setup steps for the
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
void DiffusionFluxNLFormIntegrator<VARS>::init() {
  for (const auto& u : this->u_old_) {
    this->sloth_u_old_.emplace_back(SlothGridFunction(u));
  }
  this->get_parameters();
  this->check_coefficient_types(this->expected_list_);
  this->get_coefficients();
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
void DiffusionFluxNLFormIntegrator<VARS>::AssembleElementVector(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun, const mfem::Array<mfem::Vector*>& elvect) {
  std::vector<double> vaux_gf_at_ip;
  vaux_gf_at_ip.resize(vaux_gf_.size());
  int num_blocks = el.Size();
  std::vector<double> u_values(2 * num_blocks);

  int blk = 0;
  int nd = el[blk]->GetDof();
  int dim = el[blk]->GetDim();
  gradPsi.SetSize(nd, dim);
  Psi.SetSize(nd);
  int nElement = Tr.ElementNo;

  this->Flux_.SetSize(dim);

  elvect[blk]->SetSize(nd);
  *elvect[blk] = 0.;

  mfem::Vector grad_uold;
  grad_uold.SetSize(dim);

  const mfem::IntegrationRule* ir =
      &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());

  // Loop over integration points
  for (int i = 0; i < ir->GetNPoints(); i++) {
    const mfem::IntegrationPoint& ip = ir->IntPoint(i);

    el[blk]->CalcShape(ip, Psi);
    Tr.SetIntPoint(&ip);

    // Get aux values at ip TODO(cci) (move in method)
    vaux_gf_at_ip.clear();
    for (const auto& aux_gf : vaux_gf_) {
      vaux_gf_at_ip.emplace_back(std::move(aux_gf.GetValue(Tr, ip)));
    }
    // Get values
    for (int off_blk = 0; off_blk < num_blocks; ++off_blk) {
      u_values[off_blk] = (*elfun[off_blk]) * Psi;
      u_values[off_blk + num_blocks] = this->u_old_[off_blk].GetValue(Tr, ip);
    }

    // Stabilization contribution : D_stab * (Grad u - Grad un)
    el[blk]->CalcPhysDShape(Tr, this->gradPsi);
    this->gradPsi.MultTranspose(*elfun[blk], this->Flux_);
    this->sloth_u_old_[blk].GetGradient(Tr, this->gradPsi, grad_uold);

    this->Flux_.Add(-1, grad_uold);
    this->Flux_ *= this->compute_coefficient(stab_diffusion[blk], std::span<const double>(u_values),
                                             std::span<const double>(vaux_gf_at_ip));

    // Diffusion flux (see child classes)
    this->add_diffusion_flux(Tr, nElement, ip, dim);
    double weight_coef = ip.weight * Tr.Weight();
    if (this->isAxisymmetric()) {
      weight_coef *= mfem::CylindricalRadialCoefficient().Eval(Tr, ip);
    }
    this->Flux_ *= weight_coef;

    this->gradPsi.AddMult(this->Flux_, *elvect[blk], 1.0);
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
void DiffusionFluxNLFormIntegrator<VARS>::AssembleElementGrad(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun, const mfem::Array2D<mfem::DenseMatrix*>& elmat) {
  std::vector<double> vaux_gf_at_ip;
  vaux_gf_at_ip.resize(vaux_gf_.size());
  int num_blocks = el.Size();
  std::vector<double> u_values(2 * num_blocks);
  int blk = 0;
  int nd = el[blk]->GetDof();
  int dim = el[blk]->GetDim();

  Psi.SetSize(nd);
  gradPsi.SetSize(nd, dim);

  elmat(blk, blk)->SetSize(nd);
  *elmat(blk, blk) = 0.0;
  mfem::Vector vec;
  vec.SetSize(nd);
  vec = 0.0;

  const mfem::IntegrationRule* ir =
      &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());

  for (int i = 0; i < ir->GetNPoints(); i++) {
    const mfem::IntegrationPoint& ip = ir->IntPoint(i);
    el[blk]->CalcShape(ip, Psi);
    // Get aux values at ip TODO(cci) (move in method)
    vaux_gf_at_ip.clear();
    for (const auto& aux_gf : vaux_gf_) {
      vaux_gf_at_ip.emplace_back(std::move(aux_gf.GetValue(Tr, ip)));
    }
    // Get values
    for (int off_blk = 0; off_blk < num_blocks; ++off_blk) {
      u_values[off_blk] = (*elfun[off_blk]) * Psi;
      u_values[off_blk + num_blocks] = this->u_old_[off_blk].GetValue(Tr, ip);
    }

    Tr.SetIntPoint(&ip);

    double weight_coef = ip.weight * Tr.Weight();
    if (this->isAxisymmetric()) {
      weight_coef *= mfem::CylindricalRadialCoefficient().Eval(Tr, ip);
    }

    const double coeff_diffu =
        this->compute_coefficient(stab_diffusion[blk], std::span<const double>(u_values),
                                  std::span<const double>(vaux_gf_at_ip)) *
        weight_coef;
    el[blk]->CalcPhysDShape(Tr, gradPsi);
    AddMult_a_AAt(coeff_diffu, gradPsi, *elmat(blk, blk));
  }
}

/**
 * @brief Compute the diffusion flux at the integration point
 *
 * @remark the coefficients and the gradients contributions must be diffusion in the child classes
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param Tr Element transformation.
 * @param nElement number of the element
 * @param ip Integration point
 * @param dim Spatial dimension
 */
template <class VARS>
void DiffusionFluxNLFormIntegrator<VARS>::add_diffusion_flux(mfem::ElementTransformation& Tr,
                                                             const int nElement,
                                                             const mfem::IntegrationPoint& ip,
                                                             const int dim) {
  std::vector<mfem::Vector> gradient = this->get_flux_gradient(Tr, nElement, ip, dim);
  std::vector<double> coef = this->get_flux_coefficient(nElement, ip);
  for (unsigned int i = 0; i < gradient.size(); i++) {
    this->Flux_.Add(coef[i], gradient[i]);
  }
}

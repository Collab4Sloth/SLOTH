/**
 * @file DiffusionNLFormIntegrator.tpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Nonlinear form integrator for diffusion-type problems
 * @version 0.1
 * @date 2025-09-05
 *
 * Copyright CEA (C) 2025
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
#include <vector>

#include "Coefficients/SlothBaseCoefficient.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new DiffusionNLFormIntegrator object.
 *
 * This constructor initializes a nonlinear form integrator for diffusion-type
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
 * @note This integrator serves as a base class for specific diffusion
 *       integrators like FickNLFormIntegrator and FourierNLFormIntegrator.
 */
template <class VARS>
DiffusionNLFormIntegrator<VARS>::DiffusionNLFormIntegrator(
    const std::vector<mfem::ParGridFunction> u_old,
    const std::vector<mfem::ParGridFunction> aux_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients)
    : SlothNLFormIntegrator<VARS>(u_old, aux_old, params, auxvars, coefficients) {}

/**
 * @brief Initialize the diffusion integrator.
 *
 * This method performs all necessary setup steps for a diffusion-type
 * nonlinear form integrator:
 * 1. Verifies that the list of expected coefficient types is not empty.
 * 2. Checks that all coefficients contain the expected types.
 * 3. Retrieves and stores the diffusion coefficients internally.
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 *
 * @pre The integrator's 'expected_list_' must be populated with the
 *      required coefficient types for this diffusion integrator.
 *
 * @post Internal diffusion coefficient storage is initialized.
 *
 * @throws mfem::Error If the expected list of coefficients is empty.
 */
template <class VARS>
void DiffusionNLFormIntegrator<VARS>::init() {
  MFEM_VERIFY(
      !this->expected_list_.empty(),
      "Expected not empty list of coefficients for diffusion integrators. Please check your data.");
  this->check_coefficient_types(this->expected_list_);
  this->get_coefficients();

  for (std::size_t i = 0; i < this->aux_infos_.size(); ++i) {
    this->vaux_gf_.emplace_back(std::move(this->aux_gf_[i]));
  }
}

/**
 * @brief Assemble the element-level residual vector for the nonlinear problem.
 *
 * This method computes the residual vector of the diffusion-type nonlinear problem by element
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
void DiffusionNLFormIntegrator<VARS>::AssembleElementVector(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun, const mfem::Array<mfem::Vector*>& elvect) {
  int num_blocks = el.Size();
  std::vector<double> u_values(2 * num_blocks);

  std::vector<double> vaux_gf_at_ip;
  vaux_gf_at_ip.resize(vaux_gf_.size());

  for (int blk = 0; blk < num_blocks; ++blk) {
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
      el[blk]->CalcShape(ip, Psi);
      Tr.SetIntPoint(&ip);

      // Get aux values at ip TODO(cci) (move in method)
      for (size_t k = 0; k < vaux_gf_.size(); ++k) {
        vaux_gf_at_ip[k] = vaux_gf_[k].GetValue(Tr, ip);
      }
      // Get values
      for (int off_blk = 0; off_blk < num_blocks; ++off_blk) {
        u_values[off_blk] = (*elfun[off_blk]) * Psi;
        u_values[off_blk + num_blocks] = this->u_old_[off_blk].GetValue(Tr, ip);
      }

      el[blk]->CalcPhysDShape(Tr, gradPsi);
      gradPsi.MultTranspose(*elfun[blk], gradU);
      double diffu = this->compute_coefficient(diffusion[blk], std::span<const double>(u_values),
                                               std::span<const double>(vaux_gf_at_ip));
      const double coeff_diffu = diffu * ip.weight * Tr.Weight();
      gradU *= coeff_diffu;
      gradPsi.AddMult(gradU, *elvect[blk]);
    }
  }
}
/**
 * @brief Assemble the element-level Jacobian matrix for the nonlinear problem.
 *
 * This method computes the Jacobian matrix of a diffusion-type nonlinear problem by element.
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
void DiffusionNLFormIntegrator<VARS>::AssembleElementGrad(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun, const mfem::Array2D<mfem::DenseMatrix*>& elmat) {
  int num_blocks = el.Size();
  std::vector<double> u_values(2 * num_blocks);

  std::vector<double> vaux_gf_at_ip;
  vaux_gf_at_ip.resize(vaux_gf_.size());

  for (int blk = 0; blk < num_blocks; ++blk) {
    int nd = el[blk]->GetDof();
    int dim = el[blk]->GetDim();

    Psi.SetSize(nd);
    gradPsi.SetSize(nd, dim);

    mfem::Vector vec;
    vec.SetSize(nd);

    elmat(blk, blk)->SetSize(nd);
    *elmat(blk, blk) = 0.0;

    const mfem::IntegrationRule* ir =
        &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());

    vec = 0.0;
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      el[blk]->CalcShape(ip, Psi);

      Tr.SetIntPoint(&ip);

      // Get aux values at ip TODO(cci) (move in method)
      for (size_t k = 0; k < vaux_gf_.size(); ++k) {
        vaux_gf_at_ip[k] = vaux_gf_[k].GetValue(Tr, ip);
      }
      // Get values
      for (int off_blk = 0; off_blk < num_blocks; ++off_blk) {
        u_values[off_blk] = (*elfun[off_blk]) * Psi;
        u_values[off_blk + num_blocks] = this->u_old_[off_blk].GetValue(Tr, ip);
      }

      double diffu = this->compute_coefficient(diffusion[blk], std::span<const double>(u_values),
                                               std::span<const double>(vaux_gf_at_ip));
      double grad_diffu =
          this->compute_gradient_coefficient(diffusion[blk], blk, std::span<const double>(u_values),
                                             std::span<const double>(vaux_gf_at_ip));

      el[blk]->CalcPhysDShape(Tr, gradPsi);
      AddMult_a_AAt(diffu * ip.weight * Tr.Weight(), gradPsi, *elmat(blk, blk));

      gradPsi.MultTranspose(*elfun[blk], gradU);
      gradPsi.AddMult(gradU, vec);
      AddMult_a_VWt(grad_diffu * ip.weight * Tr.Weight(), Psi, vec, *elmat(blk, blk));
    }
  }
}

/**
 * @brief Compute the value of a diffusion coefficient.
 *
 * This method evaluates the given coefficient using the provided values.
 * By default, the 'values' vector contains {u, u_old}, and auxiliary
 * variables remain accessible via the class member 'aux_gf_'.
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param coef   Coefficient to evaluate.
 * @param values Vector of current and previous solution values (default: {u, u_old}).
 *
 * @return The computed scalar value of the diffusion coefficient.
 *
 * @note This function can be overridden in derived classes.
 */
template <class VARS>
double DiffusionNLFormIntegrator<VARS>::compute_coefficient(
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
 * @brief Compute the value of a specific component of the gradient of a diffusion coefficient.
 *
 * This method evaluates the gradient of the given coefficient with respect to the
 * variable corresponding to the specified block. By default, the 'values' vector
 * contains {u, u_old}, and auxiliary variables remain accessible via the class member 'aux_gf_'.
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param coef   Coefficient whose gradient is to be computed.
 * @param iblk    Index of the gradient.
 * @param values Vector of current and previous solution values (default: {u, u_old}).
 *
 * @return The computed scalar value of the gradient component of the diffusion coefficient.
 *
 * @note This function can be overridden in derived classes.
 */
template <class VARS>
double DiffusionNLFormIntegrator<VARS>::compute_gradient_coefficient(
    Coefficient coef, const int blk, const std::span<const double>& values,
    const std::span<const double>& aux_values) {
  if (coef.is_scalar()) {
    return 0.0;
  } else {
    std::span<const double> u(values.begin(), values.begin() + this->nb_blk_);
    std::span<const double> un(values.begin() + this->nb_blk_, values.end());

    const std::span<const double>& input = coef.is_implicit() ? u : un;
    return coef.compute_gradient(blk, input, aux_values);
  }
}

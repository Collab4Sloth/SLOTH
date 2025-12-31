/**
 * @file DiffusionNLFormIntegrator.hpp
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
#include <algorithm>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "Integrators/SlothNLFormIntegrator.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

/**
 * @brief  Nonlinear form integrator for diffusion-type problems
 *
 * @tparam DIFFU_NAME
 */
template <class VARS>
class DiffusionNLFormIntegrator : public SlothNLFormIntegrator<VARS> {
 private:
  mfem::DenseMatrix gradPsi;
  mfem::Vector Psi, gradU;

 protected:
  std::list<GlossaryType> expected_list_;
  Coefficients diffusion;
  virtual double compute_coefficient(Coefficient coef, const std::vector<double>& values);
  virtual double compute_gradient_coefficient(Coefficient coef, const int blk,
                                              const std::vector<double>& values);
  void get_coefficients() override = 0;
  void init() override;

 public:
  DiffusionNLFormIntegrator(const std::vector<mfem::ParGridFunction> u_old,
                            const Parameters& params, std::vector<VARS*> auxvars,
                            const std::vector<Coefficients>& coefficients);

  void AssembleElementVector(const mfem::Array<const mfem::FiniteElement*>& el,
                             mfem::ElementTransformation& Tr,
                             const mfem::Array<const mfem::Vector*>& elfun,
                             const mfem::Array<mfem::Vector*>& elvect) override;

  void AssembleElementGrad(const mfem::Array<const mfem::FiniteElement*>& el,
                           mfem::ElementTransformation& Tr,
                           const mfem::Array<const mfem::Vector*>& elfun,
                           const mfem::Array2D<mfem::DenseMatrix*>& elmat) override;
};
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

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
    const std::vector<mfem::ParGridFunction> u_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients)
    : SlothNLFormIntegrator<VARS>(u_old, params, auxvars, coefficients) {}

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
  int blk = 0;
  int nd = el[blk]->GetDof();
  int dim = el[blk]->GetDim();
  gradPsi.SetSize(nd, dim);
  Psi.SetSize(nd);
  gradU.SetSize(dim);

  elvect[blk]->SetSize(nd);
  *elvect[blk] = 0.;

  Coefficients coeff_blk = this->coefficients_[blk];

  const mfem::IntegrationRule* ir =
      &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());

  for (int i = 0; i < ir->GetNPoints(); i++) {
    const mfem::IntegrationPoint& ip = ir->IntPoint(i);
    el[blk]->CalcShape(ip, Psi);
    Tr.SetIntPoint(&ip);

    const auto& u = *elfun[blk] * Psi;
    const auto& un = this->u_old_[blk].GetValue(Tr, ip);

    el[blk]->CalcPhysDShape(Tr, gradPsi);
    gradPsi.MultTranspose(*elfun[blk], gradU);
    double diffu = this->compute_coefficient(diffusion[blk], {u, un});
    const double coeff_diffu = diffu * ip.weight * Tr.Weight();
    gradU *= coeff_diffu;
    gradPsi.AddMult(gradU, *elvect[blk]);
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
  int blk = 0;
  int nd = el[blk]->GetDof();
  int dim = el[blk]->GetDim();

  Psi.SetSize(nd);
  gradPsi.SetSize(nd, dim);

  mfem::Vector vec;
  vec.SetSize(nd);

  elmat(blk, blk)->SetSize(nd);
  *elmat(blk, blk) = 0.0;

  Coefficients coeff_blk = this->coefficients_[blk];

  const mfem::IntegrationRule* ir =
      &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());

  vec = 0.0;
  for (int i = 0; i < ir->GetNPoints(); i++) {
    const mfem::IntegrationPoint& ip = ir->IntPoint(i);
    el[blk]->CalcShape(ip, Psi);
    const auto& u = *elfun[blk] * Psi;

    Tr.SetIntPoint(&ip);

    const auto& un = this->u_old_[blk].GetValue(Tr, ip);
    double diffu = this->compute_coefficient(diffusion[blk], {u, un});
    double grad_diffu = this->compute_gradient_coefficient(diffusion[blk], blk, {u});

    el[blk]->CalcPhysDShape(Tr, gradPsi);
    AddMult_a_AAt(diffu * ip.weight * Tr.Weight(), gradPsi, *elmat(blk, blk));

    gradPsi.MultTranspose(*elfun[blk], gradU);
    gradPsi.AddMult(gradU, vec);
    AddMult_a_VWt(grad_diffu * ip.weight * Tr.Weight(), Psi, vec, *elmat(blk, blk));
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
double DiffusionNLFormIntegrator<VARS>::compute_coefficient(Coefficient coef,
                                                            const std::vector<double>& values) {
  const double u = values[0];
  const double un = values[1];
  double coef_value = 0.0;
  if (coef.is_implicit()) {
    coef_value = coef.compute({u});
  } else if (coef.is_explicit()) {
    coef_value = coef.compute({un});
  } else if (coef.is_scalar()) {
    coef_value = coef.compute();
  }
  return coef_value;
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
    Coefficient coef, const int iblk, const std::vector<double>& values) {
  const double u = values[0];
  double coef_value = 0.0;
  if (coef.is_implicit()) {
    coef_value = coef.compute_gradient(iblk, {u});
  }
  return coef_value;
}

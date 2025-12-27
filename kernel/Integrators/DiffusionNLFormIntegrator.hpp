/**
 * @file DiffusionNLFormIntegrator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief FV for a diffusion equation
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
 * @brief  Class dedicated to the FV of the mass diffusion equation
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
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Construct a new DiffusionNLFormIntegrator<VARS>::DiffusionNLFormIntegratorobject
 *
 * @tparam VARS
 * @param u_old
 * @param params
 * @param auxvars
 * @param coefficients
 */
template <class VARS>
DiffusionNLFormIntegrator<VARS>::DiffusionNLFormIntegrator(
    const std::vector<mfem::ParGridFunction> u_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients)
    : SlothNLFormIntegrator<VARS>(u_old, params, auxvars, coefficients) {}

template <class VARS>
void DiffusionNLFormIntegrator<VARS>::init() {
  MFEM_VERIFY(
      !this->expected_list_.empty(),
      "Expected not empty list of coefficients for diffusion integrators. Please check your data.");
  this->check_coefficient_types(this->expected_list_);
  this->get_coefficients();
}

/**
 * @brief Residual part of the non linear problem
 *
 * @tparam VARS
 * @param el
 * @param Tr
 * @param elfun
 * @param elvect
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

    // Laplacian : given u, compute (grad(u), grad(psi)), psi is shape function.
    // given u (elfun), compute grad(u)
    el[blk]->CalcPhysDShape(Tr, gradPsi);
    gradPsi.MultTranspose(*elfun[blk], gradU);
    double diffu = this->compute_coefficient(diffusion[blk], {u, un});
    const double coeff_diffu = diffu * ip.weight * Tr.Weight();
    gradU *= coeff_diffu;
    gradPsi.AddMult(gradU, *elvect[blk]);
  }
}

/**
 * @brief Jacobian part of the non linear problem
 *
 * @tparam VARS
 * @param el
 * @param Tr
 * @param elfun
 * @param elmat
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
    // Laplacian : compute (grad(phi), grad(psi)), phi is shape function.

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
 * @brief Return the value of the diffusion coefficient
 * @remark by default values = {u,un} and aux_variables remain accessible in the method with the
 * class variable aux_gf_
 * @tparam VARS
 * @param coef
 * @param values
 * @return double
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
 * @brief Return the value of the component blk of the gradient of the diffusion coefficient
 * @remark by default values = {u,un} and aux_variables remain accessible in the method with the
 * class variable aux_gf_
 *
 * @tparam VARS
 * @param coef
 * @param blk
 * @param values
 * @return double
 */
template <class VARS>
double DiffusionNLFormIntegrator<VARS>::compute_gradient_coefficient(
    Coefficient coef, const int blk, const std::vector<double>& values) {
  const double u = values[0];
  const double un = values[1];
  double coef_value = 0.0;
  if (coef.is_implicit()) {
    coef_value = coef.compute_gradient(blk, {u});
  }
  return coef_value;
}

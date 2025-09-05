/**
 * @file DiffusionNLFormIntegrator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief FV for the mass diffusion equation
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

#include "Coefficients/DiffusionCoefficient.hpp"
#include "Coefficients/PhaseFieldPotentials.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

/**
 * @brief  Class dedicated to the FV of the mass diffusion equation
 *
 * @tparam SCHEME
 * @tparam DIFFU_NAME
 */
template <class VARS, CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME>
class DiffusionNLFormIntegrator : public mfem::BlockNonlinearFormIntegrator,
                                  public SlothNLFormIntegrator<VARS> {
 private:
  std::vector<mfem::ParGridFunction> u_old_;
  std::vector<mfem::ParGridFunction> aux_gf_;
  std::vector<mfem::Vector> aux_old_gf_;
  std::vector<std::vector<std::string>> aux_infos_;
  mfem::DenseMatrix gradPsi;
  mfem::Vector Psi, gradU;

  template <typename... Args>
  double diffusion(mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip,
                   const double u, const Parameters& parameters);
  template <typename... Args>
  double diffusion_prime(mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip,
                         const double u, const Parameters& parameters);

 public:
  DiffusionNLFormIntegrator(const std::vector<mfem::ParGridFunction>& u_old,
                            const Parameters& params, std::vector<VARS*> auxvars);
  ~DiffusionNLFormIntegrator();

  virtual void AssembleElementVector(const mfem::Array<const mfem::FiniteElement*>& el,
                                     mfem::ElementTransformation& Tr,
                                     const mfem::Array<const mfem::Vector*>& elfun,
                                     const mfem::Array<mfem::Vector*>& elvect);

  virtual void AssembleElementGrad(const mfem::Array<const mfem::FiniteElement*>& el,
                                   mfem::ElementTransformation& Tr,
                                   const mfem::Array<const mfem::Vector*>& elfun,
                                   const mfem::Array2D<mfem::DenseMatrix*>& elmat);

  std::unique_ptr<HomogeneousEnergyCoefficient<ThermodynamicsPotentials::LOG>> get_energy(
      std::vector<mfem::ParGridFunction*> gfu, const double diffu);
};
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Return the value of the diffusion coefficient at integration point
 *
 * @tparam SCHEME
 * @tparam DIFFU_NAME
 * @tparam Args
 * @param Tr
 * @param ip
 * @param gfu
 * @param parameters
 * @return double
 */
template <class VARS, CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME>
template <typename... Args>
double DiffusionNLFormIntegrator<VARS, SCHEME, DIFFU_NAME>::diffusion(
    mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip, const double u,
    const Parameters& parameters) {
  if (SCHEME == CoefficientDiscretization::Implicit) {
    DiffusionCoefficient<0, DIFFU_NAME> diff_coeff(u, parameters);
    return diff_coeff.Eval(Tr, ip);
  } else {
    DiffusionCoefficient<0, DIFFU_NAME> diff_coeff(&this->u_old_[0], parameters);
    return diff_coeff.Eval(Tr, ip);
  }
}

/**
 * @brief Return the first derivative of the diffusion coefficient vs the variable at integration
 * point
 *
 * @tparam SCHEME
 * @tparam DIFFU_NAME
 * @tparam Args
 * @param Tr
 * @param ip
 * @param gfu
 * @param parameters
 * @return double
 */
template <class VARS, CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME>
template <typename... Args>
double DiffusionNLFormIntegrator<VARS, SCHEME, DIFFU_NAME>::diffusion_prime(
    mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip, const double u,
    const Parameters& parameters) {
  double coef = 0.;
  if (SCHEME == CoefficientDiscretization::Implicit) {
    DiffusionCoefficient<1, DIFFU_NAME> diff_coeff(u, parameters);
    coef = diff_coeff.Eval(Tr, ip);
  }
  return coef;
}

/**
 * @brief Construct a new DiffusionNLFormIntegrator<SCHEME, COEFFICIENT>::DiffusionNLFormIntegrator
 * object
 *
 * @tparam SCHEME
 * @tparam COEFFICIENT
 * @param u_old
 * @param alpha
 * @param kappa
 */
template <class VARS, CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME>
DiffusionNLFormIntegrator<VARS, SCHEME, DIFFU_NAME>::DiffusionNLFormIntegrator(
    const std::vector<mfem::ParGridFunction>& u_old, const Parameters& params,
    std::vector<VARS*> auxvars)
    : SlothNLFormIntegrator<VARS>(params, auxvars), u_old_(u_old) {
  this->aux_gf_ = this->get_aux_gf();
  this->aux_old_gf_ = this->get_aux_old_gf();
  this->aux_infos_ = this->get_aux_infos();
}

/**
 * @brief Residual part of the non linear problem
 *
 * @param el
 * @param Tr
 * @param elfun
 * @param elvect
 */

template <class VARS, CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME>
void DiffusionNLFormIntegrator<VARS, SCHEME, DIFFU_NAME>::AssembleElementVector(
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

  const mfem::IntegrationRule* ir =
      &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());

  for (int i = 0; i < ir->GetNPoints(); i++) {
    const mfem::IntegrationPoint& ip = ir->IntPoint(i);
    el[blk]->CalcShape(ip, Psi);
    Tr.SetIntPoint(&ip);

    const auto& u = *elfun[blk] * Psi;

    // Laplacian : given u, compute (grad(u), grad(psi)), psi is shape function.
    // given u (elfun), compute grad(u)
    el[blk]->CalcPhysDShape(Tr, gradPsi);
    gradPsi.MultTranspose(*elfun[blk], gradU);

    const double coeff_diffu = this->diffusion(Tr, ip, u, this->params_) * ip.weight * Tr.Weight();
    gradU *= coeff_diffu;
    gradPsi.AddMult(gradU, *elvect[blk]);
  }
}

/**
 * @brief  Jacobian part of the non linear problem
 *
 * @param el
 * @param Tr
 * @param elfun
 * @param elmat
 */
template <class VARS, CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME>
void DiffusionNLFormIntegrator<VARS, SCHEME, DIFFU_NAME>::AssembleElementGrad(
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

  const mfem::IntegrationRule* ir =
      &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());


  vec = 0.0;
  for (int i = 0; i < ir->GetNPoints(); i++) {
    const mfem::IntegrationPoint& ip = ir->IntPoint(i);
    el[blk]->CalcShape(ip, Psi);
    const auto& u = *elfun[blk] * Psi;
    // Laplacian : compute (grad(phi), grad(psi)), phi is shape function.

    Tr.SetIntPoint(&ip);
    const double coeff_diffu = this->diffusion(Tr, ip, u, this->params_) * ip.weight * Tr.Weight();
    el[blk]->CalcPhysDShape(Tr, gradPsi);
    AddMult_a_AAt(coeff_diffu, gradPsi, *elmat(blk, blk));

    gradPsi.MultTranspose(*elfun[blk], gradU);
    gradPsi.AddMult(gradU, vec);
    const auto coef_diffu_derivative = this->diffusion_prime(Tr, ip, u, this->params_);
    AddMult_a_VWt(coef_diffu_derivative * ip.weight * Tr.Weight(), Psi, vec, *elmat(blk, blk));
  }
}

template <class VARS, CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME>
std::unique_ptr<HomogeneousEnergyCoefficient<ThermodynamicsPotentials::LOG>>
DiffusionNLFormIntegrator<VARS, SCHEME, DIFFU_NAME>::get_energy(
    std::vector<mfem::ParGridFunction*> gfu, const double diffu) {
  return std::make_unique<HomogeneousEnergyCoefficient<ThermodynamicsPotentials::LOG>>(gfu[0],
                                                                                       diffu);
}

/**
 * @brief Destroy the DiffusionNLFormIntegrator<SCHEME, COEFFICIENT>::DiffusionNLFormIntegrator
 *
 * @tparam SCHEME
 * @tparam COEFFICIENT
 */

template <class VARS, CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME>
DiffusionNLFormIntegrator<VARS, SCHEME, DIFFU_NAME>::~DiffusionNLFormIntegrator() {}

/**
 * @file DiffusionFluxNLFormIntegrator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief inter-diffusion integrator
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

#include <memory>
#include <utility>
#include <vector>

#include "Integrators/SlothGridFunction.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameters.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

/**
 * @brief  Class dedicated to the VF of an inter-diffusion equation
 *
 * @tparam SCHEME
 * @tparam DIFFU_NAME
 */
template <class VARS>
class DiffusionFluxNLFormIntegrator : public mfem::BlockNonlinearFormIntegrator,
                                      public SlothNLFormIntegrator<VARS> {
 private:
  double coeff_stab_;

  void add_diffusion_flux(mfem::ElementTransformation& Tr, const int nElement,
                          const mfem::IntegrationPoint& ip, const int dim);

 protected:
  std::vector<SlothGridFunction> u_old_;

  mfem::DenseMatrix gradPsi;
  mfem::Vector Psi, Flux_;

  virtual void get_parameters();

  virtual std::vector<mfem::Vector> get_flux_gradient(mfem::ElementTransformation& Tr,
                                                      const int nElement,
                                                      const mfem::IntegrationPoint& ip,
                                                      const int dim) = 0;
  virtual std::vector<double> get_flux_coefficient(const int nElement,
                                                   const mfem::IntegrationPoint& ip) = 0;

 public:
  DiffusionFluxNLFormIntegrator(const std::vector<mfem::ParGridFunction>& u_old,
                                const Parameters& params, std::vector<VARS*> auxvars);
  ~DiffusionFluxNLFormIntegrator();

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
 * @brief Return parameters required by these integrators
 *
 * @tparam VARS
 */
template <class VARS>
void DiffusionFluxNLFormIntegrator<VARS>::get_parameters() {
  // Get stabilization coefficient
  this->coeff_stab_ = this->params_.template get_param_value<double>("D");
}

/**
 * @brief Construct a new DiffusionFluxNLFormIntegrator<VARS>::DiffusionFluxNLFormIntegrator object
 *
 * @tparam VARS
 * @param u_old
 * @param params
 * @param auxvars
 */
template <class VARS>
DiffusionFluxNLFormIntegrator<VARS>::DiffusionFluxNLFormIntegrator(
    const std::vector<mfem::ParGridFunction>& u_old, const Parameters& params,
    std::vector<VARS*> auxvars)
    : SlothNLFormIntegrator<VARS>(params, auxvars) {
  for (const auto& u : u_old) {
    this->u_old_.emplace_back(std::move(SlothGridFunction(u)));
  }
  this->get_parameters();
}

/**
 * @brief Residual part
 *
 * @tparam VARS
 * @param el
 * @param Tr
 * @param elfun
 * @param elvect
 */
template <class VARS>
void DiffusionFluxNLFormIntegrator<VARS>::AssembleElementVector(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun, const mfem::Array<mfem::Vector*>& elvect) {
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

    const auto& u = *elfun[blk] * Psi;
    // Stabilization contribution : D_stab * (Grad u - Grad un)
    el[blk]->CalcPhysDShape(Tr, this->gradPsi);
    this->gradPsi.MultTranspose(*elfun[blk], this->Flux_);
    this->u_old_[blk].GetGradient(Tr, this->gradPsi, grad_uold);

    this->Flux_.Add(-1, grad_uold);
    this->Flux_ *= this->coeff_stab_;

    // Diffusion flux (see child classes)
    this->add_diffusion_flux(Tr, nElement, ip, dim);

    this->Flux_ *= ip.weight * Tr.Weight();

    this->gradPsi.AddMult(this->Flux_, *elvect[blk], 1.0);
  }
}

/**
 * @brief Jacobian part
 *
 * @tparam VARS
 * @param el
 * @param Tr
 * @param elfun
 * @param elmat
 */
template <class VARS>
void DiffusionFluxNLFormIntegrator<VARS>::AssembleElementGrad(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun, const mfem::Array2D<mfem::DenseMatrix*>& elmat) {
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
    const auto& u = *elfun[blk] * Psi;

    Tr.SetIntPoint(&ip);
    const double coeff_diffu = this->coeff_stab_ * ip.weight * Tr.Weight();
    el[blk]->CalcPhysDShape(Tr, gradPsi);
    AddMult_a_AAt(coeff_diffu, gradPsi, *elmat(blk, blk));
  }
}

/**
 * @brief Return the diffusion flux
 * @remark the coefficients and the gradients contributions must be diffusion in the child classes
 *
 * @tparam VARS
 * @param Tr
 * @param nElement
 * @param ip
 * @param dim
 */
template <class VARS>
void DiffusionFluxNLFormIntegrator<VARS>::add_diffusion_flux(mfem::ElementTransformation& Tr,
                                                             const int nElement,
                                                             const mfem::IntegrationPoint& ip,
                                                             const int dim) {
  std::vector<mfem::Vector> gradient = this->get_flux_gradient(Tr, nElement, ip, dim);
  std::vector<double> coef = this->get_flux_coefficient(nElement, ip);
  for (int i = 0; i < gradient.size(); i++) {
    this->Flux_.Add(coef[i], gradient[i]);
  }
}

/**
 * @brief
 *
 * @tparam VARS
 * @param gfu
 * @param diffu
 * @return std::unique_ptr<HomogeneousEnergyCoefficient<ThermodynamicsPotentials::LOG>>
 */
template <class VARS>
std::unique_ptr<HomogeneousEnergyCoefficient<ThermodynamicsPotentials::LOG>>
DiffusionFluxNLFormIntegrator<VARS>::get_energy(std::vector<mfem::ParGridFunction*> gfu,
                                                const double diffu) {
  return std::make_unique<HomogeneousEnergyCoefficient<ThermodynamicsPotentials::LOG>>(gfu[0],
                                                                                       diffu);
}

/**
 * @brief Destroy the DiffusionFluxNLFormIntegrator< VARS>::DiffusionFluxNLFormIntegrator object
 *
 * @tparam VARS
 */
template <class VARS>
DiffusionFluxNLFormIntegrator<VARS>::~DiffusionFluxNLFormIntegrator() {}

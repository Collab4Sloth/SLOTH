/**
 * @file DiffusionFluxNLFormIntegrator.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief inter-diffusion integrator
 * @version 0.1
 * @date 2025-04-03
 *
 * @copyright Copyright (c) 2025
 *
 */

#include <algorithm>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Coefficients/DiffusionCoefficient.hpp"
#include "Coefficients/PhaseFieldPotentials.hpp"
#include "Integrators/SlothGridFunction.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
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
class DiffusionFluxNLFormIntegrator : public mfem::NonlinearFormIntegrator,
                                      public SlothNLFormIntegrator<VARS> {
 private:
  double coeff_stab_;

  void add_diffusion_flux(mfem::ElementTransformation& Tr, const int nElement,
                          const mfem::IntegrationPoint& ip, const int dim);

 protected:
  SlothGridFunction u_old_;

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
  DiffusionFluxNLFormIntegrator(const mfem::ParGridFunction& u_old, const Parameters& params,
                                std::vector<VARS*> auxvars);
  ~DiffusionFluxNLFormIntegrator();

  virtual void AssembleElementVector(const mfem::FiniteElement& el, mfem::ElementTransformation& Tr,
                                     const mfem::Vector& elfun, mfem::Vector& elvect);

  virtual void AssembleElementGrad(const mfem::FiniteElement& el, mfem::ElementTransformation& Tr,
                                   const mfem::Vector& elfun, mfem::DenseMatrix& elmat);

  std::unique_ptr<HomogeneousEnergyCoefficient<ThermodynamicsPotentials::LOG>> get_energy(
      mfem::ParGridFunction* gfu, const double diffu);
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
    const mfem::ParGridFunction& u_old, const Parameters& params, std::vector<VARS*> auxvars)
    : SlothNLFormIntegrator<VARS>(params, auxvars), u_old_(u_old) {
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
void DiffusionFluxNLFormIntegrator<VARS>::AssembleElementVector(const mfem::FiniteElement& el,
                                                                mfem::ElementTransformation& Tr,
                                                                const mfem::Vector& elfun,
                                                                mfem::Vector& elvect) {
  int nd = el.GetDof();
  int dim = el.GetDim();
  int nElement = Tr.ElementNo;

  this->Psi.SetSize(nd);
  this->gradPsi.SetSize(nd, dim);

  this->Flux_.SetSize(dim);

  elvect.SetSize(nd);
  mfem::Vector grad_uold;
  grad_uold.SetSize(dim);

  // Initialization
  elvect = 0.0;

  const mfem::IntegrationRule* ir =
      &mfem::IntRules.Get(el.GetGeomType(), 2 * el.GetOrder() + Tr.OrderW());

  // Loop over integration points
  for (int i = 0; i < ir->GetNPoints(); i++) {
    const mfem::IntegrationPoint& ip = ir->IntPoint(i);

    el.CalcShape(ip, Psi);
    Tr.SetIntPoint(&ip);

    const auto& u = elfun * Psi;

    // Stabilization contribution : D_stab * (Grad u - Grad un)
    el.CalcPhysDShape(Tr, this->gradPsi);
    this->gradPsi.MultTranspose(elfun, this->Flux_);
    this->u_old_.GetGradient(Tr, this->gradPsi, grad_uold);

    this->Flux_.Add(-1, grad_uold);
    this->Flux_ *= this->coeff_stab_;

    // Diffusion flux (see child classes)
    this->add_diffusion_flux(Tr, nElement, ip, dim);

    this->Flux_ *= ip.weight * Tr.Weight();

    this->gradPsi.AddMult(this->Flux_, elvect, 1.0);
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
void DiffusionFluxNLFormIntegrator<VARS>::AssembleElementGrad(const mfem::FiniteElement& el,
                                                              mfem::ElementTransformation& Tr,
                                                              const mfem::Vector& elfun,
                                                              mfem::DenseMatrix& elmat) {
  int nd = el.GetDof();
  int dim = el.GetDim();

  Psi.SetSize(nd);
  gradPsi.SetSize(nd, dim);
  elmat.SetSize(nd);
  mfem::Vector vec;
  vec.SetSize(nd);

  const mfem::IntegrationRule* ir =
      &mfem::IntRules.Get(el.GetGeomType(), 2 * el.GetOrder() + Tr.OrderW());

  elmat = 0.0;
  vec = 0.0;
  for (int i = 0; i < ir->GetNPoints(); i++) {
    const mfem::IntegrationPoint& ip = ir->IntPoint(i);
    el.CalcShape(ip, Psi);
    const auto& u = elfun * Psi;

    Tr.SetIntPoint(&ip);
    const double coeff_diffu = this->coeff_stab_ * ip.weight * Tr.Weight();
    el.CalcPhysDShape(Tr, gradPsi);
    AddMult_a_AAt(coeff_diffu, gradPsi, elmat);
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
  std::vector<double> coef = this->get_flux_coefficient(Tr, nElement, ip);
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
DiffusionFluxNLFormIntegrator<VARS>::get_energy(mfem::ParGridFunction* gfu, const double diffu) {
  return std::make_unique<HomogeneousEnergyCoefficient<ThermodynamicsPotentials::LOG>>(gfu, diffu);
}

/**
 * @brief Destroy the DiffusionFluxNLFormIntegrator< VARS>::DiffusionFluxNLFormIntegrator object
 *
 * @tparam VARS
 */
template <class VARS>
DiffusionFluxNLFormIntegrator<VARS>::~DiffusionFluxNLFormIntegrator() {}

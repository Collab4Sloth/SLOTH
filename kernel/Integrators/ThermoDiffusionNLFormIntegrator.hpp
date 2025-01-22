/**
 * @file ThermoDiffusionNLFormIntegrator.hpp
 * @author cp273896  (clement.plumecocq@cea.fr)
 * @brief VF of an inter-diffusion equation
 * @version 0.1
 * @date 2024-06-06
 *
 * Copyright CEA (c) 2024
 *
 */
#include <algorithm>
#include <memory>
#include <string>
#include <tuple>
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
template <class VARS, CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME>
class ThermoDiffusionNLFormIntegrator : public mfem::NonlinearFormIntegrator,
                                        public SlothNLFormIntegrator<VARS> {
 private:
  std::vector<std::tuple<std::string, double>> inter_diffusion_coeff_;
  double coeff_stab_;

  SlothGridFunction u_old_;

  std::vector<mfem::ParGridFunction> mu_gf_;
  mfem::DenseMatrix gradPsi;
  mfem::Vector Psi, gradMu_;

 protected:
  void get_parameters(const Parameters& vectr_param);

 public:
  ThermoDiffusionNLFormIntegrator(const mfem::ParGridFunction& u_old, const Parameters& params,
                                  std::vector<VARS*> auxvars);
  ~ThermoDiffusionNLFormIntegrator();

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
template <class VARS, CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME>
void ThermoDiffusionNLFormIntegrator<VARS, SCHEME, DIFFU_NAME>::get_parameters(
    const Parameters& params) {
  MFEM_VERIFY(
      this->mu_gf_.size() == (params.get_size() - 1),
      "ThermoDiffusionNLFormIntegrator requires as many parameters (inter-diffusion "
      "coefficients) as auxiliary variables in addition of the stabilization coefficient D");

  this->coeff_stab_ = params.get_param_value<double>("D");

  for (const auto& p : params.get_vector()) {
    const std::string& para_name = p.get_name();
    if (para_name == "D") continue;
    const auto& value = p.get_value();
    const double& para_value = std::get<double>(value);
    this->inter_diffusion_coeff_.emplace_back(para_name, para_value);
  }
}

/**
 * @brief Construct a new ThermoDiffusionNLFormIntegrator<SCHEME,
 * COEFFICIENT>::ThermoDiffusionNLFormIntegrator object
 *
 * @tparam SCHEME
 * @tparam COEFFICIENT
 * @param u_old
 * @param alpha
 * @param kappa
 */
template <class VARS, CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME>
ThermoDiffusionNLFormIntegrator<VARS, SCHEME, DIFFU_NAME>::ThermoDiffusionNLFormIntegrator(
    const mfem::ParGridFunction& u_old, const Parameters& params, std::vector<VARS*> auxvars)
    : SlothNLFormIntegrator<VARS>(params, auxvars), u_old_(u_old) {
  this->mu_gf_ = this->get_auxiliary_gf();
  this->get_parameters(params);
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
void ThermoDiffusionNLFormIntegrator<VARS, SCHEME, DIFFU_NAME>::AssembleElementVector(
    const mfem::FiniteElement& el, mfem::ElementTransformation& Tr, const mfem::Vector& elfun,
    mfem::Vector& elvect) {
  int nd = el.GetDof();
  int dim = el.GetDim();
  int nElement = Tr.ElementNo;

  this->Psi.SetSize(nd);
  this->gradPsi.SetSize(nd, dim);

  this->gradMu_.SetSize(dim);

  elvect.SetSize(nd);
  mfem::Vector grad_uold;
  mfem::Vector grad_mu_i;
  grad_mu_i.SetSize(dim);
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
    this->gradPsi.MultTranspose(elfun, this->gradMu_);
    this->u_old_.GetGradient(Tr, this->gradPsi, grad_uold);

    this->gradMu_ -= grad_uold;
    this->gradMu_ *= this->coeff_stab_;

    // Sum of [M_i * Grad mu_i] : explicit contribution (not present in Jacobian)
    for (size_t i = 0; i < this->mu_gf_.size(); i++) {
      auto mu_i = SlothGridFunction(this->mu_gf_.at(i));
      mu_i.GetGradient(Tr, this->gradPsi, grad_mu_i);

      const auto& M_i = std::get<1>(this->inter_diffusion_coeff_[i]) *
                        this->u_old_.GetValue(nElement, ip) *
                        (1. - this->u_old_.GetValue(nElement, ip));
      grad_mu_i *= M_i;
      this->gradMu_ += grad_mu_i;
    }
    this->gradMu_ *= ip.weight * Tr.Weight();

    this->gradPsi.AddMult(this->gradMu_, elvect, 1.0);
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
void ThermoDiffusionNLFormIntegrator<VARS, SCHEME, DIFFU_NAME>::AssembleElementGrad(
    const mfem::FiniteElement& el, mfem::ElementTransformation& Tr, const mfem::Vector& elfun,
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

template <class VARS, CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME>
std::unique_ptr<HomogeneousEnergyCoefficient<ThermodynamicsPotentials::LOG>>
ThermoDiffusionNLFormIntegrator<VARS, SCHEME, DIFFU_NAME>::get_energy(mfem::ParGridFunction* gfu,
                                                                      const double diffu) {
  return std::make_unique<HomogeneousEnergyCoefficient<ThermodynamicsPotentials::LOG>>(gfu, diffu);
}

/**
 * @brief Destroy the ThermoDiffusionNLFormIntegrator<SCHEME,
 * COEFFICIENT>::ThermoDiffusionNLFormIntegrator
 *
 * @tparam SCHEME
 * @tparam COEFFICIENT
 */

template <class VARS, CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME>
ThermoDiffusionNLFormIntegrator<VARS, SCHEME, DIFFU_NAME>::~ThermoDiffusionNLFormIntegrator() {}

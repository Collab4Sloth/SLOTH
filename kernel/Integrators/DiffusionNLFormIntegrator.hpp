/**
 * @file DiffusionNLFormIntegrator.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief FV for the mass diffusion equation
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
class DiffusionNLFormIntegrator : public mfem::NonlinearFormIntegrator,
                                  public SlothNLFormIntegrator<VARS> {
 private:
  mfem::ParGridFunction u_old_;
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
  DiffusionNLFormIntegrator(const mfem::ParGridFunction& u_old, const Parameters& params,
                            std::vector<VARS*> auxvars);
  ~DiffusionNLFormIntegrator();

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
    DiffusionCoefficient<0, DIFFU_NAME> diff_coeff(&this->u_old_, parameters);
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
    const mfem::ParGridFunction& u_old, const Parameters& params, std::vector<VARS*> auxvars)
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
    const mfem::FiniteElement& el, mfem::ElementTransformation& Tr, const mfem::Vector& elfun,
    mfem::Vector& elvect) {
  int nd = el.GetDof();
  int dim = el.GetDim();
  gradPsi.SetSize(nd, dim);
  Psi.SetSize(nd);
  gradU.SetSize(dim);
  elvect.SetSize(nd);
  const mfem::IntegrationRule* ir =
      &mfem::IntRules.Get(el.GetGeomType(), 2 * el.GetOrder() + Tr.OrderW());
  elvect = 0.0;
  for (int i = 0; i < ir->GetNPoints(); i++) {
    const mfem::IntegrationPoint& ip = ir->IntPoint(i);
    el.CalcShape(ip, Psi);
    Tr.SetIntPoint(&ip);

    const auto& u = elfun * Psi;

    // Laplacian : given u, compute (grad(u), grad(psi)), psi is shape function.
    // given u (elfun), compute grad(u)
    el.CalcPhysDShape(Tr, gradPsi);
    gradPsi.MultTranspose(elfun, gradU);

    const double coeff_diffu = this->diffusion(Tr, ip, u, this->params_) * ip.weight * Tr.Weight();
    gradU *= coeff_diffu;
    gradPsi.AddMult(gradU, elvect);
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
    // Laplacian : compute (grad(phi), grad(psi)), phi is shape function.

    Tr.SetIntPoint(&ip);
    const double coeff_diffu = this->diffusion(Tr, ip, u, this->params_) * ip.weight * Tr.Weight();
    el.CalcPhysDShape(Tr, gradPsi);
    AddMult_a_AAt(coeff_diffu, gradPsi, elmat);

    gradPsi.MultTranspose(elfun, gradU);
    gradPsi.AddMult(gradU, vec);
    const auto coef_diffu_derivative = this->diffusion_prime(Tr, ip, u, this->params_);
    AddMult_a_VWt(coef_diffu_derivative * ip.weight * Tr.Weight(), Psi, vec, elmat);
  }
}

template <class VARS, CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME>
std::unique_ptr<HomogeneousEnergyCoefficient<ThermodynamicsPotentials::LOG>>
DiffusionNLFormIntegrator<VARS, SCHEME, DIFFU_NAME>::get_energy(mfem::ParGridFunction* gfu,
                                                                const double diffu) {
  return std::make_unique<HomogeneousEnergyCoefficient<ThermodynamicsPotentials::LOG>>(gfu, diffu);
}

/**
 * @brief Destroy the DiffusionNLFormIntegrator<SCHEME, COEFFICIENT>::DiffusionNLFormIntegrator
 *
 * @tparam SCHEME
 * @tparam COEFFICIENT
 */

template <class VARS, CoefficientDiscretization SCHEME, Diffusion DIFFU_NAME>
DiffusionNLFormIntegrator<VARS, SCHEME, DIFFU_NAME>::~DiffusionNLFormIntegrator() {}

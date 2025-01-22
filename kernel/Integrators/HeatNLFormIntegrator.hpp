/**
 * @file HeatNLFormIntegrator.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief FV for the heat transfer equation
 * @version 0.1
 * @date 2024-06-06
 *
 * Copyright CEA (c) 2024
 *
 */
#include <algorithm>
#include <memory>
#include <tuple>
#include <vector>

#include "Coefficients/ConductivityCoefficient.hpp"
#include "Coefficients/PhaseFieldPotentials.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

/**
 * @brief Class dedicated to the FV of the heat transfer equation
 *
 * @tparam SCHEME
 * @tparam COND_NAME
 */
template <class VARS, CoefficientDiscretization SCHEME, Conductivity COND_NAME>
class HeatNLFormIntegrator : public mfem::NonlinearFormIntegrator,
                             public SlothNLFormIntegrator<VARS> {
 private:
  mfem::ParGridFunction u_old_;
  std::vector<mfem::ParGridFunction> aux_gf_;
  mfem::DenseMatrix gradPsi;
  mfem::Vector Psi, gradU;

  template <typename... Args>
  double conductivity(mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip,
                      const double u, const Parameters& parameters);
  template <typename... Args>
  double conductivity_prime(mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip,
                            const double u, const Parameters& parameters);

 public:
  HeatNLFormIntegrator(const mfem::ParGridFunction& u_old, const Parameters& params,
                       std::vector<VARS*> auxvars);
  ~HeatNLFormIntegrator();

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
 * @brief Return the value of the Conductivity coefficient at integration point
 *
 * @tparam SCHEME
 * @tparam COND_NAME
 * @tparam Args
 * @param Tr
 * @param ip
 * @param gfu
 * @param parameters
 * @return double
 */
template <class VARS, CoefficientDiscretization SCHEME, Conductivity COND_NAME>
template <typename... Args>
double HeatNLFormIntegrator<VARS, SCHEME, COND_NAME>::conductivity(mfem::ElementTransformation& Tr,
                                                                   const mfem::IntegrationPoint& ip,
                                                                   const double u,
                                                                   const Parameters& parameters) {
  if (SCHEME == CoefficientDiscretization::Implicit) {
    ConductivityCoefficient<0, COND_NAME> coeff(u, parameters);
    return coeff.Eval(Tr, ip);
  } else {
    ConductivityCoefficient<0, COND_NAME> coeff(&this->u_old_, parameters);
    return coeff.Eval(Tr, ip);
  }
}

/**
 * @brief Return the first derivative of the Conductivity coefficient vs the variable at integration
 * point
 *
 * @tparam SCHEME
 * @tparam COND_NAME
 * @tparam Args
 * @param Tr
 * @param ip
 * @param gfu
 * @param parameters
 * @return double
 */
template <class VARS, CoefficientDiscretization SCHEME, Conductivity COND_NAME>
template <typename... Args>
double HeatNLFormIntegrator<VARS, SCHEME, COND_NAME>::conductivity_prime(
    mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip, const double u,
    const Parameters& parameters) {
  double coef = 0.;
  if (SCHEME == CoefficientDiscretization::Implicit) {
    ConductivityCoefficient<1, COND_NAME> cond_coeff(u, parameters);
    coef = cond_coeff.Eval(Tr, ip);
  }
  return coef;
}

/**
 * @brief Construct a new HeatNLFormIntegrator<SCHEME, COEFFICIENT>::HeatNLFormIntegrator
 * object
 *
 * @tparam SCHEME
 * @tparam COEFFICIENT
 * @param u_old
 * @param alpha
 * @param kappa
 */
template <class VARS, CoefficientDiscretization SCHEME, Conductivity COND_NAME>
HeatNLFormIntegrator<VARS, SCHEME, COND_NAME>::HeatNLFormIntegrator(
    const mfem::ParGridFunction& u_old, const Parameters& params, std::vector<VARS*> auxvars)
    : SlothNLFormIntegrator<VARS>(params, auxvars), u_old_(u_old) {
  this->aux_gf_ = this->get_aux_gf();
}

/**
 * @brief Residual part of the non linear problem
 *
 * @param el
 * @param Tr
 * @param elfun
 * @param elvect
 */

template <class VARS, CoefficientDiscretization SCHEME, Conductivity COND_NAME>
void HeatNLFormIntegrator<VARS, SCHEME, COND_NAME>::AssembleElementVector(
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

    const double coeff_diffu =
        this->conductivity(Tr, ip, u, this->params_) * ip.weight * Tr.Weight();
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
template <class VARS, CoefficientDiscretization SCHEME, Conductivity COND_NAME>
void HeatNLFormIntegrator<VARS, SCHEME, COND_NAME>::AssembleElementGrad(
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
    const double coeff_diffu =
        this->conductivity(Tr, ip, u, this->params_) * ip.weight * Tr.Weight();
    el.CalcPhysDShape(Tr, gradPsi);
    AddMult_a_AAt(coeff_diffu, gradPsi, elmat);

    gradPsi.MultTranspose(elfun, gradU);
    gradPsi.AddMult(gradU, vec);
    const auto coef_diffu_derivative = this->conductivity_prime(Tr, ip, u, this->params_);
    AddMult_a_VWt(coef_diffu_derivative * ip.weight * Tr.Weight(), Psi, vec, elmat);
  }
}

template <class VARS, CoefficientDiscretization SCHEME, Conductivity COND_NAME>
std::unique_ptr<HomogeneousEnergyCoefficient<ThermodynamicsPotentials::LOG>>
HeatNLFormIntegrator<VARS, SCHEME, COND_NAME>::get_energy(mfem::ParGridFunction* gfu,
                                                          const double diffu) {
  return std::make_unique<HomogeneousEnergyCoefficient<ThermodynamicsPotentials::LOG>>(gfu, diffu);
}

/**
 * @brief Destroy the HeatNLFormIntegrator<SCHEME, COEFFICIENT>::HeatNLFormIntegrator
 *
 * @tparam SCHEME
 * @tparam COEFFICIENT
 */

template <class VARS, CoefficientDiscretization SCHEME, Conductivity COND_NAME>
HeatNLFormIntegrator<VARS, SCHEME, COND_NAME>::~HeatNLFormIntegrator() {}

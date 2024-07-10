/**
 * @file DiffusionNLFormIntegrator.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief
 * @version 0.1
 * @date 2024-06-06
 *
 * Copyright CEA (c) 2024
 *
 */
#include <algorithm>
#include <tuple>

#include "Coefficients/DiffusionCoeffients.hpp"
#include "Coefficients/SourceTermCoefficient.hpp"
#include "Utils/PhaseFieldOptions.hpp"
#include "mfem.hpp" // NOLINT [no include the directory when naming mfem include file]

#pragma once

using FuncType = std::function<double(const double&, const double&)>;

template <DiffusionCoefficientDiscretization SCHEME, DiffusionCoefficients COEFFICIENT>
class DiffusionNLFormIntegrator : public mfem::NonlinearFormIntegrator {
 private:
  mfem::GridFunction u_old_;
  mfem::DenseMatrix gradPsi;
  mfem::Vector Psi, gradU;

  DiffusionFunctions<0, SCHEME, COEFFICIENT> diffusion_coefficient_function_;
  DiffusionFunctions<1, SCHEME, COEFFICIENT> diffusion_coefficient_function_first_derivative_;

  template <typename... Args>
  FuncType diffusion(const int order_derivative, Args... parameters);

 protected:
  double alpha_, kappa_;

 public:
  DiffusionNLFormIntegrator(const mfem::GridFunction& u_old, const double& alpha,
                            const double& kappa);
  ~DiffusionNLFormIntegrator();

  virtual void AssembleElementVector(const mfem::FiniteElement& el, mfem::ElementTransformation& Tr,
                                     const mfem::Vector& elfun, mfem::Vector& elvect);

  virtual void AssembleElementGrad(const mfem::FiniteElement& el, mfem::ElementTransformation& Tr,
                                   const mfem::Vector& elfun, mfem::DenseMatrix& elmat);
};
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Diffusion coefficient (or derivative)
 *
 * @tparam SCHEME
 * @tparam COEFFICIENT
 * @param order_derivative
 * @return FuncType
 */
template <DiffusionCoefficientDiscretization SCHEME, DiffusionCoefficients COEFFICIENT>
template <typename... Args>
FuncType DiffusionNLFormIntegrator<SCHEME, COEFFICIENT>::diffusion(const int order_derivative,
                                                                   Args... parameters) {
  return FuncType([this, order_derivative, parameters...](const double& u, const double& un) {
    std::function<double(const double&)> func;
    if (order_derivative == 0) {
      func = this->diffusion_coefficient_function_.getFunction(un, parameters...);
    } else if (order_derivative == 1) {
      func = this->diffusion_coefficient_function_first_derivative_.getFunction(un, parameters...);
    } else {
      throw std::runtime_error(
          "Error while setting the order of derivative: only 0 and 1 are allowed.");
    }

    return func(u);
  });
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
template <DiffusionCoefficientDiscretization SCHEME, DiffusionCoefficients COEFFICIENT>
DiffusionNLFormIntegrator<SCHEME, COEFFICIENT>::DiffusionNLFormIntegrator(
    const mfem::GridFunction& u_old, const double& alpha, const double& kappa)
    : u_old_(u_old), alpha_(alpha), kappa_(kappa) {}

/**
 * @brief Residual part of the non linear problem
 *
 * @param el
 * @param Tr
 * @param elfun
 * @param elvect
 */
template <DiffusionCoefficientDiscretization SCHEME, DiffusionCoefficients COEFFICIENT>
void DiffusionNLFormIntegrator<SCHEME, COEFFICIENT>::AssembleElementVector(
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
    const auto& un = this->u_old_.GetValue(Tr, ip);

    // Laplacian : given u, compute (grad(u), grad(psi)), psi is shape function.
    // given u (elfun), compute grad(u)
    el.CalcPhysDShape(Tr, gradPsi);
    gradPsi.MultTranspose(elfun, gradU);

    const double coeff_diffu =
        this->diffusion(0, this->alpha_, this->kappa_)(u, un) * ip.weight * Tr.Weight();
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
template <DiffusionCoefficientDiscretization SCHEME, DiffusionCoefficients COEFFICIENT>
void DiffusionNLFormIntegrator<SCHEME, COEFFICIENT>::AssembleElementGrad(
    const mfem::FiniteElement& el, mfem::ElementTransformation& Tr, const mfem::Vector& elfun,
    mfem::DenseMatrix& elmat) {
  int nd = el.GetDof();
  int dim = el.GetDim();
  int spaceDim = Tr.GetSpaceDim();
  bool square = (dim == spaceDim);
  double w;

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
    const auto& un = this->u_old_.GetValue(Tr, ip);
    // Laplacian : compute (grad(phi), grad(psi)), phi is shape function.
    Tr.SetIntPoint(&ip);
    const double coeff_diffu =
        this->diffusion(0, this->alpha_, this->kappa_)(u, un) * ip.weight * Tr.Weight();
    el.CalcPhysDShape(Tr, gradPsi);
    AddMult_a_AAt(coeff_diffu, gradPsi, elmat);

    gradPsi.MultTranspose(elfun, gradU);
    gradPsi.AddMult(gradU, vec);
    const auto coef_diffu_derivative = this->diffusion(1, this->alpha_, this->kappa_)(u, un);
    AddMult_a_VWt(coef_diffu_derivative * ip.weight * Tr.Weight(), Psi, vec, elmat);
  }
}

/**
 * @brief Destroy the DiffusionNLFormIntegrator<SCHEME, COEFFICIENT>::DiffusionNLFormIntegrator
 *
 * @tparam SCHEME
 * @tparam COEFFICIENT
 */
template <DiffusionCoefficientDiscretization SCHEME, DiffusionCoefficients COEFFICIENT>
DiffusionNLFormIntegrator<SCHEME, COEFFICIENT>::~DiffusionNLFormIntegrator() {}

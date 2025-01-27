/**
 * @file AllenCahnNLFormIntegrator.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief FV for the Allen Cahn equation
 * @todo Extend coefficient to omega and lambda
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

#include "Coefficients/LambdaCoefficient.hpp"
#include "Coefficients/MobilityCoefficient.hpp"
#include "Coefficients/OmegaCoefficient.hpp"
#include "Coefficients/PhaseFieldPotentials.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Profiling/Profiling.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

/**
 * @brief Class dedicated to the FV of the Allen Cahn equation
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
class AllenCahnNLFormIntegrator : public mfem::NonlinearFormIntegrator,
                                  public SlothNLFormIntegrator<VARS> {
 private:
  mfem::DenseMatrix gradPsi;
  mfem::Vector Psi, gradU;

  PotentialFunctions<1, SCHEME, ENERGY> energy_first_derivative_potential_;
  PotentialFunctions<2, SCHEME, ENERGY> energy_second_derivative_potential_;

  template <typename... Args>
  double mobility(mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip, const double u,
                  const Parameters& parameters);

  template <typename... Args>
  double omega(mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip, const double u,
               const Parameters& parameters);

  template <typename... Args>
  double lambda(mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip, const double u,
                const Parameters& parameters);

  FType double_well_derivative(const int order_derivative, mfem::ElementTransformation& Tr,
                               const mfem::IntegrationPoint& ir);

 protected:
  mfem::ParGridFunction u_old_;
  std::vector<mfem::ParGridFunction> aux_gf_;

  virtual FType energy_derivatives(const int order_derivative, mfem::ElementTransformation& Tr,
                                   const mfem::IntegrationPoint& ir);

 public:
  AllenCahnNLFormIntegrator(const mfem::ParGridFunction& u_old, const Parameters& params,
                            std::vector<VARS*> auxvars);
  ~AllenCahnNLFormIntegrator();

  virtual void AssembleElementVector(const mfem::FiniteElement& el, mfem::ElementTransformation& Tr,
                                     const mfem::Vector& elfun, mfem::Vector& elvect);

  virtual void AssembleElementGrad(const mfem::FiniteElement& el, mfem::ElementTransformation& Tr,
                                   const mfem::Vector& elfun, mfem::DenseMatrix& elmat);

  std::unique_ptr<HomogeneousEnergyCoefficient<ENERGY>> get_energy(mfem::ParGridFunction* gfu,
                                                                   const double omega);
  std::unique_ptr<GradientEnergyCoefficient> get_grad_energy(mfem::ParGridFunction* gfu,
                                                             const double lambda);
};
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Energy derivatives contribution
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @param order_derivative
 * @return FuncType
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
FType AllenCahnNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::energy_derivatives(
    const int order_derivative, mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ir) {
  return [this, order_derivative, &Tr, &ir](const double& u) {
    return this->double_well_derivative(order_derivative, Tr, ir)(u);
  };
}

/**
 * @brief Return the value of the mobility coefficient at integration point
 *
 * @remark Actually, only explicit mobility is available. (See diffusion integrator  to extend
 * towards implicit mobility. )
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam Args
 * @param Tr
 * @param ip
 * @param u
 * @param parameters
 * @return double
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
template <typename... Args>
double AllenCahnNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::mobility(
    mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip, const double u,
    const Parameters& parameters) {
  MobilityCoefficient<0, MOBI> mobi_coeff(&this->u_old_, parameters);
  return mobi_coeff.Eval(Tr, ip);
}

/**
 * @brief Return the value of the omega coefficient at integration point
 *
 * @remark Actually, only constant omega is available.
 *
 * @tparam VARS
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam Args
 * @param Tr
 * @param ip
 * @param u
 * @param parameters
 * @return double
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
template <typename... Args>
double AllenCahnNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::omega(
    mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip, const double u,
    const Parameters& parameters) {
  OmegaCoefficient<0, Omega::Constant> omega_coeff(&this->u_old_, parameters);
  return omega_coeff.Eval(Tr, ip);
}
/**
 * @brief Return the value of the lambda coefficient at integration point
 *
 * @remark Actually, only constant lambda is available.
 *
 * @tparam VARS
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam Args
 * @param Tr
 * @param ip
 * @param u
 * @param parameters
 * @return double
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
template <typename... Args>
double AllenCahnNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::lambda(
    mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip, const double u,
    const Parameters& parameters) {
  LambdaCoefficient<0, Lambda::Constant> lambda_coeff(&this->u_old_, parameters);
  return lambda_coeff.Eval(Tr, ip);
}

/**
 * @brief double_well_derivative contribution
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @param order_derivative
 * @return std::function<double(const double&, const double&)>
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
FType AllenCahnNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::double_well_derivative(
    const int order_derivative, mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ir) {
  return FType([this, order_derivative, &Tr, &ir](const double& u) {
    const auto& un = this->u_old_.GetValue(Tr, ir);
    FType W_derivative;
    if (order_derivative == 1) {
      W_derivative = this->energy_first_derivative_potential_.getPotentialFunction(un);
    } else if (order_derivative == 2) {
      W_derivative = this->energy_second_derivative_potential_.getPotentialFunction(un);
    } else {
      std::runtime_error("Error while setting the order of derivative : only 1 and 2 are allowed.");
    }
    const double omega = this->omega(Tr, ir, u, this->params_);

    const auto& w_prime = omega * W_derivative(u);
    return w_prime;
  });
}

/**
 * @brief Construct a new AllenCahnNLFormIntegrator object
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @param u_old
 * @param omega
 * @param lambda
 * @param mob
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
AllenCahnNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::AllenCahnNLFormIntegrator(
    const mfem::ParGridFunction& u_old, const Parameters& params, std::vector<VARS*> auxvars)
    : SlothNLFormIntegrator<VARS>(params, auxvars), u_old_(u_old) {
  this->aux_gf_ = this->get_aux_gf();
}

/**
 * @brief Residual part of the non linear problem
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @param el
 * @param Tr
 * @param elfun
 * @param elvect
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
void AllenCahnNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::AssembleElementVector(
    const mfem::FiniteElement& el, mfem::ElementTransformation& Tr, const mfem::Vector& elfun,
    mfem::Vector& elvect) {
  Catch_Time_Section("AllenCahnNLFormIntegrator:AssembleElementVector");
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
    el.CalcShape(ip, Psi);  //
    Tr.SetIntPoint(&ip);

    const auto& u = elfun * Psi;

    // Laplacian : given u, compute (grad(u), grad(psi)), psi is shape function.
    // given u (elfun), compute grad(u)
    el.CalcPhysDShape(Tr, gradPsi);
    gradPsi.MultTranspose(elfun, gradU);
    const double coef_mob = this->mobility(Tr, ip, u, this->params_) * ip.weight * Tr.Weight();
    const double lambda = this->lambda(Tr, ip, u, this->params_);
    gradU *= coef_mob * lambda;
    gradPsi.AddMult(gradU, elvect);

    // Given u, compute (w'(u), psi), psi is shape function
    const double ww = coef_mob * this->energy_derivatives(1, Tr, ip)(u);
    add(elvect, ww, Psi, elvect);
  }
}

/**
 * @brief Jacobian part of the non linear problem
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @param el
 * @param Tr
 * @param elfun
 * @param elmat
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
void AllenCahnNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::AssembleElementGrad(
    const mfem::FiniteElement& el, mfem::ElementTransformation& Tr, const mfem::Vector& elfun,
    mfem::DenseMatrix& elmat) {
  Catch_Time_Section("AllenCahnNLFormIntegrator::AssembleElementGrad");

  int nd = el.GetDof();
  int dim = el.GetDim();

  gradPsi.SetSize(nd, dim);
  Psi.SetSize(nd);
  elmat.SetSize(nd);

  const mfem::IntegrationRule* ir =
      &mfem::IntRules.Get(el.GetGeomType(), 2 * el.GetOrder() + Tr.OrderW());

  elmat = 0.0;
  for (int i = 0; i < ir->GetNPoints(); i++) {
    const mfem::IntegrationPoint& ip = ir->IntPoint(i);
    el.CalcShape(ip, Psi);
    const auto& u = elfun * Psi;
    Tr.SetIntPoint(&ip);

    // Laplacian : compute (grad(u), grad(psi)), psi is shape function.
    const double coef_mob = this->mobility(Tr, ip, u, this->params_) * ip.weight * Tr.Weight();
    el.CalcPhysDShape(Tr, gradPsi);
    const double lambda = this->lambda(Tr, ip, u, this->params_);

    AddMult_a_AAt(coef_mob * lambda, gradPsi, elmat);

    // Compute w'(u)*(du,psi), psi is shape function ( // w''(u))
    double fun_val = coef_mob * this->energy_derivatives(2, Tr, ip)(u);
    AddMult_a_VVt(fun_val, Psi, elmat);  // w'(u)*(du, psi)
  }
}

/**
 * @brief Return the energy coefficient associated with the integrator
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam Args
 * @param gfu
 * @param lambda
 * @param omega
 * @return EnergyCoefficient<SCHEME, ENERGY>
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
std::unique_ptr<HomogeneousEnergyCoefficient<ENERGY>>
AllenCahnNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::get_energy(mfem::ParGridFunction* gfu,
                                                                  const double omega) {
  return std::make_unique<HomogeneousEnergyCoefficient<ENERGY>>(gfu, omega);
}

/**
 * @brief  Return the gradient energy coefficient associated with the integrator
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @param gfu
 * @param lambda
 * @return GradientEnergyCoefficient
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
std::unique_ptr<GradientEnergyCoefficient>
AllenCahnNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::get_grad_energy(mfem::ParGridFunction* gfu,
                                                                       const double lambda) {
  return std::make_unique<GradientEnergyCoefficient>(gfu, lambda);
}

/**
 * @brief Destroy the AllenCahnNLFormIntegrator  object
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
AllenCahnNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::~AllenCahnNLFormIntegrator() {}

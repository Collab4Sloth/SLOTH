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

#include "Coefficients/MobilityCoefficient.hpp"
#include "Coefficients/PhaseFieldPotentials.hpp"
#include "Coefficients/SourceTermCoefficient.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Profiling/Profiling.hpp"
#include "Utils/PhaseFieldOptions.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

/**
 * @brief Class dedicated to the FV of the Allen Cahn equation
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 */
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
class AllenCahnNLFormIntegrator : public mfem::NonlinearFormIntegrator {
 private:
  const Parameters mobility_params_;
  mfem::DenseMatrix dshape, dshapedxt, invdfdx;
  mfem::Vector shape, vec, pointflux;

  PotentialFunctions<1, SCHEME, ENERGY> energy_first_derivative_potential_;
  PotentialFunctions<2, SCHEME, ENERGY> energy_second_derivative_potential_;

  template <typename... Args>
  double mobility(mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip, const double u,
                  const Parameters& parameters);

  FType double_well_derivative(const int order_derivative, mfem::ElementTransformation& Tr,
                               const mfem::IntegrationPoint& ir);

 protected:
  std::vector<mfem::ParGridFunction> aux_gf_;
  mfem::ParGridFunction u_old_;
  double omega_, lambda_;

  virtual FType energy_derivatives(const int order_derivative, mfem::ElementTransformation& Tr,
                                   const mfem::IntegrationPoint& ir);
  virtual void get_parameters(const Parameters& vectr_param);

 public:
  AllenCahnNLFormIntegrator(const mfem::ParGridFunction& u_old, const Parameters& params,
                            const std::vector<mfem::ParGridFunction>& aux_gf);
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
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
FType AllenCahnNLFormIntegrator<SCHEME, ENERGY, MOBI>::energy_derivatives(
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
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
template <typename... Args>
double AllenCahnNLFormIntegrator<SCHEME, ENERGY, MOBI>::mobility(mfem::ElementTransformation& Tr,
                                                                 const mfem::IntegrationPoint& ip,
                                                                 const double u,
                                                                 const Parameters& parameters) {
  MobilityCoefficient<0, MOBI> mobi_coeff(&this->u_old_, parameters);
  return mobi_coeff.Eval(Tr, ip);
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
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
FType AllenCahnNLFormIntegrator<SCHEME, ENERGY, MOBI>::double_well_derivative(
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
    const auto& w_prime = this->omega_ * W_derivative(u);
    return w_prime;
  });
}

template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
void AllenCahnNLFormIntegrator<SCHEME, ENERGY, MOBI>::get_parameters(const Parameters& params) {
  this->omega_ = params.get_param_value<double>("omega");
  this->lambda_ = params.get_param_value<double>("lambda");
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
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
AllenCahnNLFormIntegrator<SCHEME, ENERGY, MOBI>::AllenCahnNLFormIntegrator(
    const mfem::ParGridFunction& u_old, const Parameters& params,
    const std::vector<mfem::ParGridFunction>& aux_gf)
    : u_old_(u_old), mobility_params_(params), aux_gf_(aux_gf) {
  this->get_parameters(params);
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
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
void AllenCahnNLFormIntegrator<SCHEME, ENERGY, MOBI>::AssembleElementVector(
    const mfem::FiniteElement& el, mfem::ElementTransformation& Tr, const mfem::Vector& elfun,
    mfem::Vector& elvect) {
  Catch_Time_Section("AllenCahnNLFormIntegrator:AssembleElementVector");

  int nd = el.GetDof();
  int dim = el.GetDim();
  int spaceDim = Tr.GetSpaceDim();
  dshape.SetSize(nd, dim);
  shape.SetSize(nd);
  invdfdx.SetSize(dim, spaceDim);
  vec.SetSize(dim);
  pointflux.SetSize(spaceDim);

  elvect.SetSize(nd);
  const mfem::IntegrationRule* ir =
      &mfem::IntRules.Get(el.GetGeomType(), 2 * el.GetOrder() + Tr.OrderW());
  elvect = 0.0;
  for (int i = 0; i < ir->GetNPoints(); i++) {
    const mfem::IntegrationPoint& ip = ir->IntPoint(i);
    el.CalcDShape(ip, dshape);  // dphi
    el.CalcShape(ip, shape);    // phi
    Tr.SetIntPoint(&ip);

    const auto& u = elfun * shape;
    // const auto& un = this->u_old_.GetValue(Tr, ip);

    // Given phi, compute (w'(phi), v), v is shape function
    const double& ww = this->mobility(Tr, ip, u, this->mobility_params_) * ip.weight * Tr.Weight() *
                       this->energy_derivatives(1, Tr, ip)(u);
    add(elvect, ww, shape, elvect);

    // Laplacian : given u, compute (grad(u), grad(v)), v is shape function.
    CalcAdjugate(Tr.Jacobian(), invdfdx);  // invdfdx = adj(J)
    dshape.MultTranspose(elfun, vec);
    invdfdx.MultTranspose(vec, pointflux);
    double w =
        this->mobility(Tr, ip, u, this->mobility_params_) * this->lambda_ * ip.weight / Tr.Weight();
    pointflux *= w;
    invdfdx.Mult(pointflux, vec);

    //
    dshape.AddMult(vec, elvect);
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
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
void AllenCahnNLFormIntegrator<SCHEME, ENERGY, MOBI>::AssembleElementGrad(
    const mfem::FiniteElement& el, mfem::ElementTransformation& Tr, const mfem::Vector& elfun,
    mfem::DenseMatrix& elmat) {
  Catch_Time_Section("AllenCahnNLFormIntegrator::AssembleElementGrad");

  int nd = el.GetDof();
  int dim = el.GetDim();
  int spaceDim = Tr.GetSpaceDim();
  bool square = (dim == spaceDim);
  double w;

  shape.SetSize(nd);
  dshape.SetSize(nd, dim);
  dshapedxt.SetSize(nd, spaceDim);
  elmat.SetSize(nd);

  const mfem::IntegrationRule* ir =
      &mfem::IntRules.Get(el.GetGeomType(), 2 * el.GetOrder() + Tr.OrderW());

  elmat = 0.0;
  for (int i = 0; i < ir->GetNPoints(); i++) {
    const mfem::IntegrationPoint& ip = ir->IntPoint(i);
    el.CalcDShape(ip, dshape);  // dphi
    const auto& u = elfun * shape;
    // const auto& un = this->u_old_.GetValue(Tr, ip);

    Tr.SetIntPoint(&ip);
    w = Tr.Weight();  // det(J)
    // std::cout << " SQUARE  ? " << square << std::endl;
    w = ip.weight / (square ? w : w * w * w);
    // AdjugateJacobian = / adj(J),         if J is square
    //                    \ adj(J^t.J).J^t, otherwise

    // Tr.AdjugateJacobian() det(J)J-1

    w *= this->mobility(Tr, ip, u, this->mobility_params_) * this->lambda_;

    // dshapedxt =  det(J)J-1 dshape
    Mult(dshape, Tr.AdjugateJacobian(), dshapedxt);
    // elmat += w * dshapedxt * dshapedxt^T
    AddMult_a_AAt(w, dshapedxt, elmat);

    // Compute w'(u)*(du,v), v is shape function ( // w''(u))
    double fun_val = this->mobility(Tr, ip, u, this->mobility_params_) *
                     this->energy_derivatives(2, Tr, ip)(u) * ip.weight * Tr.Weight();
    // elmat += fun_val * shape * shape^T
    AddMult_a_VVt(fun_val, shape, elmat);  // w'(u)*(du, v)
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
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
std::unique_ptr<HomogeneousEnergyCoefficient<ENERGY>>
AllenCahnNLFormIntegrator<SCHEME, ENERGY, MOBI>::get_energy(mfem::ParGridFunction* gfu,
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
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
std::unique_ptr<GradientEnergyCoefficient>
AllenCahnNLFormIntegrator<SCHEME, ENERGY, MOBI>::get_grad_energy(mfem::ParGridFunction* gfu,
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
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
AllenCahnNLFormIntegrator<SCHEME, ENERGY, MOBI>::~AllenCahnNLFormIntegrator() {}

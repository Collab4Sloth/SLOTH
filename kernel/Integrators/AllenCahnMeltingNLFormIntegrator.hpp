/*
 * AllenCahnMeltingNLFormIntegrator.h
 *
 *  Created on: 15 may 2023
 *      Author: ci230846
 */
#include <algorithm>
#include <tuple>
#include "Coefficients/MobilityCoefficient.hpp"
#include "Coefficients/PhaseChangeFunction.hpp"
#include "Coefficients/PhaseFieldMobilities.hpp"
#include "Coefficients/PhaseFieldPotentials.hpp"
#include "Utils/PhaseFieldOptions.hpp"
#include "mfem.hpp"

#pragma once
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          ThermodynamicsPotentials INTERPOLATION, Mobility MOBI, PhaseChange PHASECHANGE>
class AllenCahnMeltingNLFormIntegrator : public mfem::NonlinearFormIntegrator {
 private:
  mfem::GridFunction u_old;
  mfem::Vector shape;
  double omega, lambda;
  double mob_;
  double alpha_;

  PhaseChangeFunction<PHASECHANGE> phase_change_function_;

  PotentialFunctions<1, SCHEME, ENERGY> energy_first_derivative_potential_;
  PotentialFunctions<1, SCHEME, INTERPOLATION> interpolation_first_derivative_potential_;
  PotentialFunctions<2, SCHEME, ENERGY> energy_second_derivative_potential_;
  PotentialFunctions<2, SCHEME, INTERPOLATION> interpolation_second_derivative_potential_;

  MobilityFunctions<MOBI> mobility_function_;

  mfem::DenseMatrix dshape, dshapedxt, invdfdx;
  mfem::Vector vec, pointflux;

 public:
  AllenCahnMeltingNLFormIntegrator(const mfem::GridFunction& _u_old, const double& _omega,
                                   const double& _lambda, const double& _mob, const double& _alpha);

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
 * @brief Construct
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam INTERPOLATION
 * @param _u_old
 * @param _omega
 * @param _lambda
 * @param _alpha
 * @param _mob
 * @return AllenCahnMeltingNLFormIntegrator<SCHEME, ENERGY, INTERPOLATION>::
 */
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          ThermodynamicsPotentials INTERPOLATION, Mobility MOBI, PhaseChange PHASECHANGE>
AllenCahnMeltingNLFormIntegrator<SCHEME, ENERGY, INTERPOLATION, MOBI, PHASECHANGE>::
    AllenCahnMeltingNLFormIntegrator(const mfem::GridFunction& _u_old, const double& _omega,
                                     const double& _lambda, const double& _mob,
                                     const double& _alpha)
    : u_old(_u_old), omega(_omega), lambda(_lambda), mob_(_mob), alpha_(_alpha) {}

/**
 * @brief Residual part of the non linear problem
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam INTERPOLATION
 * @param el
 * @param Tr
 * @param elfun
 * @param elvect
 */
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          ThermodynamicsPotentials INTERPOLATION, Mobility MOBI, PhaseChange PHASECHANGE>
void AllenCahnMeltingNLFormIntegrator<SCHEME, ENERGY, INTERPOLATION, MOBI, PHASECHANGE>::
    AssembleElementVector(const mfem::FiniteElement& el, mfem::ElementTransformation& Tr,
                          const mfem::Vector& elfun, mfem::Vector& elvect) {
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

    const auto u = elfun * shape;
    const auto un = u_old.GetValue(Tr, ip);

    const auto W = this->energy_first_derivative_potential_.getPotentialFunction(un);
    const auto H = this->interpolation_first_derivative_potential_.getPotentialFunction(un);
    const auto Wprime = W(u);
    const auto Hprime = H(u);
    const auto MOB = this->mobility_function_.getMobilityFunction(un);
    const auto Mphi = MOB(this->mob_);
    const auto PHChange = this->phase_change_function_.getPhaseChangeFunction(un);
    const auto phase_change = PHChange(this->alpha_);

    CalcAdjugate(Tr.Jacobian(), invdfdx);  // invdfdx = adj(J)

    dshape.MultTranspose(elfun, vec);
    invdfdx.MultTranspose(vec, pointflux);

    // Energy contribution
    const auto energy = this->omega * Wprime;
    // Enthalpy of melting contribution
    const auto enthalpy = phase_change * Hprime;
    // flow of phase change + source term
    const auto fun_val = Mphi * (energy + enthalpy);

    // Given phi, compute (w'(phi)-f, v), v is shape function
    const double ww = ip.weight * Tr.Weight() * fun_val;
    add(elvect, ww, shape, elvect);

    // Laplacian : given u, compute (grad(u), grad(v)), v is shape function.
    double w;
    w = Mphi * ip.weight * this->lambda / Tr.Weight();
    pointflux *= w;
    invdfdx.Mult(pointflux, vec);
    dshape.AddMult(vec, elvect);
  }
}

/**
 * @brief Jacobian part of the non linear problem
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam INTERPOLATION
 * @param el
 * @param Tr
 * @param elfun
 * @param elmat
 */
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          ThermodynamicsPotentials INTERPOLATION, Mobility MOBI, PhaseChange PHASECHANGE>
void AllenCahnMeltingNLFormIntegrator<SCHEME, ENERGY, INTERPOLATION, MOBI, PHASECHANGE>::
    AssembleElementGrad(const mfem::FiniteElement& el, mfem::ElementTransformation& Tr,
                        const mfem::Vector& elfun, mfem::DenseMatrix& elmat) {
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
    const auto u = elfun * shape;
    const auto un = u_old.GetValue(Tr, ip);

    const auto W = this->energy_second_derivative_potential_.getPotentialFunction(un);
    const auto H = this->interpolation_second_derivative_potential_.getPotentialFunction(un);
    const auto Wsecond = W(u);
    const auto Hsecond = H(u);
    const auto MOB = this->mobility_function_.getMobilityFunction(un);
    const auto Mphi = MOB(this->mob_);
    const auto PHChange = this->phase_change_function_.getPhaseChangeFunction(un);
    const auto phase_change = PHChange(this->alpha_);

    Tr.SetIntPoint(&ip);
    w = Tr.Weight();  // det(J)
    // std::cout << " SQUARE  ? " << square << std::endl;
    w = ip.weight / (square ? w : w * w * w);
    // AdjugateJacobian = / adj(J),         if J is square
    //                    \ adj(J^t.J).J^t, otherwise

    // Tr.AdjugateJacobian() det(J)J-1

    // w = w* Mphi * lambda
    w *= Mphi * this->lambda;

    // dshapedxt =  det(J)J-1 dshape
    Mult(dshape, Tr.AdjugateJacobian(), dshapedxt);
    // elmat += w * dshapedxt * dshapedxt^T
    AddMult_a_AAt(w, dshapedxt, elmat);

    //  (this->omega * secondDerivativedoubleWellPotential(elfun * shape) +
    //   this->alpha * secondDerivativeInterpolationPotential(elfun * shape)) *
    // Compute w'(u)*(du,v), v is shape function
    double fun_val =
        Mphi * (this->omega * Wsecond + phase_change * Hsecond) * ip.weight * Tr.Weight();  // w'(u)
    // elmat += fun_val * shape * shape^T
    AddMult_a_VVt(fun_val, shape, elmat);  // w'(u)*(du, v)
  }
}

/**
 * @file AllenCahnNLFormIntegrator.hpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief
 * @version 0.1
 * @date 2024-05-01
 *
 * @copyright Copyright (c) 2024
 *
 */
#include <algorithm>
#include <tuple>

#include "Coefficients/MobilityCoefficient.hpp"
#include "Coefficients/PhaseFieldMobilities.hpp"
#include "Coefficients/PhaseFieldPotentials.hpp"
#include "Coefficients/SourceTermCoefficient.hpp"
#include "Profiling/UtilsforOutput.hpp"
#include "Profiling/output.hpp"
#include "Profiling/timers.hpp"
#include "Utils/PhaseFieldOptions.hpp"
#include "mfem.hpp"

template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
class AllenCahnNLFormIntegrator : public mfem::NonlinearFormIntegrator {
 private:
  mfem::GridFunction u_old;
  mfem::Vector shape;
  double omega, lambda;
  double mob_;

  PotentialFunctions<1, SCHEME, ENERGY> energy_first_derivative_potential_;
  PotentialFunctions<2, SCHEME, ENERGY> energy_second_derivative_potential_;

  MobilityFunctions<MOBI> mobility_function_;

  mfem::DenseMatrix dshape, dshapedxt, invdfdx;
  mfem::Vector vec, pointflux;

 public:
  AllenCahnNLFormIntegrator(const mfem::GridFunction& _u_old, const double& _omega,
                            const double& _lambda, const double& _mob);

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
 * @param _u_old
 * @param _omega
 * @param _lambda
 * @param _alpha
 * @param _mob
 * @return AllenCahnNLFormIntegrator<SCHEME, ENERGY>::
 */
template <ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
AllenCahnNLFormIntegrator<SCHEME, ENERGY, MOBI>::AllenCahnNLFormIntegrator(
    const mfem::GridFunction& _u_old, const double& _omega, const double& _lambda,
    const double& _mob)
    : u_old(_u_old), omega(_omega), lambda(_lambda), mob_(_mob) {}

/**
 * @brief Residual part of the non linear problem
 *
 * @tparam SCHEME
 * @tparam ENERGY
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
  //-------------------------------
  //------- PROFILING
  //-------------------------------
  Timers timer_CalcAdjugate("CalcAdjugate");
  Timers timer_MultTranspose("MultTranspose");
  Timers timer_VMult("VMult");
  Timers timer_AddMult("AddMult");
  Timers timer_AssembleElementVector("AssembleElementVector");
  //-------------------------------
  timer_AssembleElementVector.start();

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
    const auto Wprime = W(u);
    const auto MOB = this->mobility_function_.getMobilityFunction(un);
    const auto Mphi = MOB(this->mob_);

    //--------------------------
    timer_CalcAdjugate.start();
    CalcAdjugate(Tr.Jacobian(), invdfdx);  // invdfdx = adj(J)
    timer_CalcAdjugate.stop();
    UtilsForOutput::getInstance().update_timer("CalcAdjugate", timer_CalcAdjugate);
    //--------------------------

    //--------------------------
    timer_MultTranspose.start();
    dshape.MultTranspose(elfun, vec);
    invdfdx.MultTranspose(vec, pointflux);
    timer_MultTranspose.stop();
    UtilsForOutput::getInstance().update_timer("MultTranspose", timer_MultTranspose);
    //--------------------------

    // Energy contribution + source
    const auto energy = this->omega * Wprime;
    const auto fun_val = Mphi * energy;

    // Given phi, compute (w'(phi)-f, v), v is shape function
    const double ww = ip.weight * Tr.Weight() * fun_val;
    add(elvect, ww, shape, elvect);

    // Laplacian : given u, compute (grad(u), grad(v)), v is shape function.
    double w;
    w = Mphi * ip.weight * this->lambda / Tr.Weight();
    pointflux *= w;

    //--------------------------
    timer_VMult.start();
    invdfdx.Mult(pointflux, vec);
    timer_VMult.stop();
    UtilsForOutput::getInstance().update_timer("VMult", timer_VMult);
    //--------------------------

    //--------------------------
    timer_AddMult.start();
    dshape.AddMult(vec, elvect);
    timer_AddMult.stop();
    UtilsForOutput::getInstance().update_timer("AddMult", timer_AddMult);
    //--------------------------
  }

  //--------------------------
  // save the results of profiling
  timer_AssembleElementVector.stop();
  UtilsForOutput::getInstance().update_timer("AssembleElementVector", timer_AssembleElementVector);
  //--------------------------
}

/**
 * @brief Jacobian part of the non linear problem
 *
 * @tparam SCHEME
 * @tparam ENERGY
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
  //------Start profiling-------------------------
  Timers timer_AssembleElementGrad("AssembleElementGrad");
  Timers timer_AddMult_a_VVt("AddMult_a_VVt");
  Timers timer_AddMult_a_AAt("AddMult_a_AAt");
  Timers timer_Mult("Mult");
  timer_AssembleElementGrad.start();
  //------------------------------------------------

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
    const auto Wsecond = W(u);
    const auto MOB = this->mobility_function_.getMobilityFunction(un);
    const auto Mphi = MOB(this->mob_);

    Tr.SetIntPoint(&ip);
    w = Tr.Weight();  // det(J)
    // std::cout << " SQUARE  ? " << square << std::endl;
    w = ip.weight / (square ? w : w * w * w);
    // AdjugateJacobian = / adj(J),         if J is square
    //                    \ adj(J^t.J).J^t, otherwise

    // Tr.AdjugateJacobian() det(J)J-1

    w *= Mphi * this->lambda;

    // dshapedxt =  det(J)J-1 dshape
    //--------------------------
    timer_Mult.start();
    Mult(dshape, Tr.AdjugateJacobian(), dshapedxt);
    timer_Mult.stop();
    UtilsForOutput::getInstance().update_timer("Mult", timer_Mult);
    //--------------------------

    // elmat += w * dshapedxt * dshapedxt^T
    //--------------------------
    timer_AddMult_a_AAt.start();
    AddMult_a_AAt(w, dshapedxt, elmat);
    timer_AddMult_a_AAt.stop();
    UtilsForOutput::getInstance().update_timer("AddMult_a_AAt", timer_AddMult_a_AAt);
    //--------------------------

    //  (this->omega * secondDerivativedoubleWellPotential(elfun * shape) +
    //   this->alpha * secondDerivativeInterpolationPotential(elfun * shape)) *
    // Compute w'(u)*(du,v), v is shape function
    double fun_val = (Mphi * this->omega * Wsecond) * ip.weight * Tr.Weight();  // w'(u)

    // elmat += fun_val * shape * shape^T

    //--------------------------
    timer_AddMult_a_VVt.start();
    AddMult_a_VVt(fun_val, shape, elmat);  // w'(u)*(du, v)
    timer_AddMult_a_VVt.stop();
    UtilsForOutput::getInstance().update_timer("AddMult_a_VVt", timer_AddMult_a_VVt);
    //--------------------------
  }

  // save the results of profiling
  timer_AssembleElementGrad.stop();
  UtilsForOutput::getInstance().update_timer("AssembleElementGrad", timer_AssembleElementGrad);
  //------------------------------------------------
}

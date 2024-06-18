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

#include "Coefficients/MobilityCoefficient.hpp"
#include "Coefficients/PhaseFieldMobilities.hpp"
#include "Coefficients/SourceTermCoefficient.hpp"
#include "Utils/PhaseFieldOptions.hpp"
#include "mfem.hpp"

#pragma once

using FuncType = std::function<double(const double&, const double&)>;
class DiffusionNLFormIntegrator : public mfem::NonlinearFormIntegrator {
 private:
  mfem::GridFunction u_old_;
  mfem::DenseMatrix dshape, dshapedxt, invdfdx;
  mfem::Vector shape, vec, pointflux;

  FuncType laplacian();
  FuncType stabilization();

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
 * @brief Stabilization terms
 *
 * @tparam SCHEME
 * @tparam MOBI
 * @return FuncType
 */

FuncType DiffusionNLFormIntegrator::stabilization() {
  return [](const double& u, const double& un) { return 0.; };
}

/**
 * @brief Laplacian coefficient
 *
 * @tparam SCHEME
 * @tparam MOBI
 * @return std::function<double(const double&, const double&)>
 */

FuncType DiffusionNLFormIntegrator::laplacian() {
  return FuncType([this](const double& u, const double& un) {
    const auto& laplacian = this->alpha_ * un + this->kappa_;
    return laplacian;
  });
}

/**
 * @brief Construct a new DiffusionNLFormIntegrator object
 *
 * @param u_old
 * @param alpha
 * @param kappa
 */
DiffusionNLFormIntegrator::DiffusionNLFormIntegrator(const mfem::GridFunction& u_old,
                                                     const double& alpha, const double& kappa)
    : u_old_(u_old), alpha_(alpha), kappa_(kappa) {}

/**
 * @brief Residual part of the non linear problem
 *
 * @param el
 * @param Tr
 * @param elfun
 * @param elvect
 */
void DiffusionNLFormIntegrator::AssembleElementVector(const mfem::FiniteElement& el,
                                                      mfem::ElementTransformation& Tr,
                                                      const mfem::Vector& elfun,
                                                      mfem::Vector& elvect) {
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
    const auto& un = this->u_old_.GetValue(Tr, ip);
    // mfem::Vector gradun;
    // this->u_old_.GetGradient(Tr, gradun);

    // Given phi, compute (w'(phi), v), v is shape function
    const double& ww = 0.;
    // ip.weight* Tr.Weight() * this->stabilization()(u, un);
    add(elvect, ww, shape, elvect);

    // Laplacian : given u, compute (grad(u), grad(v)), v is shape function.
    CalcAdjugate(Tr.Jacobian(), invdfdx);   // invdfdx = adj(J)
    dshape.MultTranspose(elfun, vec);       //= dphi^t*elfun = dphi^t*u
    invdfdx.MultTranspose(vec, pointflux);  //= Jadj^t*dphi^t*u : grad phi sur l'element K

    double w;
    w = this->laplacian()(u, un) * ip.weight / Tr.Weight();  // wq*D/det(J)
    pointflux *= w;                                          //= wq/det(J)*rD*Q(xq)*Jadj^t*dphi*u
    invdfdx.Mult(pointflux, vec);

    //
    dshape.AddMult(vec, elvect);
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
void DiffusionNLFormIntegrator::AssembleElementGrad(const mfem::FiniteElement& el,
                                                    mfem::ElementTransformation& Tr,
                                                    const mfem::Vector& elfun,
                                                    mfem::DenseMatrix& elmat) {
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
    const auto& un = this->u_old_.GetValue(Tr, ip);

    Tr.SetIntPoint(&ip);
    w = Tr.Weight();  // det(J)
    // std::cout << " SQUARE  ? " << square << std::endl;
    w = ip.weight / (square ? w : w * w * w);
    // AdjugateJacobian = / adj(J),         if J is square
    //                    \ adj(J^t.J).J^t, otherwise

    // Tr.AdjugateJacobian() det(J)J-1

    w *= this->laplacian()(u, un);

    // dshapedxt =  det(J)J-1 dshape
    Mult(dshape, Tr.AdjugateJacobian(), dshapedxt);
    // elmat += w * dshapedxt * dshapedxt^T
    AddMult_a_AAt(w, dshapedxt, elmat);

    // Compute w'(u)*(du,v), v is shape function ( // w''(u))
    // double fun_val = this->stabilization()(u, un) * ip.weight * Tr.Weight();
    // elmat += fun_val * shape * shape^T
    // AddMult_a_VVt(fun_val, shape, elmat);  // w'(u)*(du, v)
  }
}

/**
 * @brief Destroy the DiffusionNLFormIntegrator  object
 *
 */
DiffusionNLFormIntegrator::~DiffusionNLFormIntegrator() {}

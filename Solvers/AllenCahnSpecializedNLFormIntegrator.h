/*
 * AllenCahnSpecializedNLFormIntegrator.h
 *
 *  Created on: 8 déc. 2021
 *      Author: ci230846
 */
#include "mfem.hpp"

#ifndef SOLVERS_ALLENCAHNSPECIALIZEDNLFORMINTEGRATOR_H_
#define SOLVERS_ALLENCAHNSPECIALIZEDNLFORMINTEGRATOR_H_

class AllenCahnSpecializedNLFormIntegrator
    : public mfem::NonlinearFormIntegrator {
 private:
  mfem::Vector shape;
  mfem::Coefficient* Calphad;  // f in F(u)=-Laplace u + u^2 - q
  mfem::ConstantCoefficient *mob, *lambda, *omega;
  mfem::DenseMatrix dshape, dshapedxt, invdfdx;
  mfem::Vector vec, pointflux;

  double doubleWellPotential(const double& x) {
    return x * x * (1.0 - x) * (1.0 - x);
  };
  double firstDerivativedoubleWellPotential(const double& x) {
    return x * (1.0 - x) * (1.0 - 2 * x);
  };
  double secondDerivativedoubleWellPotential(const double& x) {
    return (1.0 - 6. * x * +6. * x * x);
  };

 public:
  AllenCahnSpecializedNLFormIntegrator(mfem::Coefficient& _Calphad,
                                       mfem::ConstantCoefficient& _mob,
                                       mfem::ConstantCoefficient& _lambda,
                                       mfem::ConstantCoefficient& _omega)
      : Calphad(&_Calphad), mob(&_mob), lambda(&_lambda), omega(&_omega){};

  virtual void AssembleElementVector(const mfem::FiniteElement& el,
                                     mfem::ElementTransformation& Tr,
                                     const mfem::Vector& elfun,
                                     mfem::Vector& elvect) {
    // Linearized DoubleWell contribution
    // W' = W'(phi_n)+W''(phi_n)*(phi_n+1 - phi_n)
    // with :
    // W = phi^2 * (1-phi)^2
    // W' = 2 phi * (1-phi) * (1-2*phi)
    // w'' = 2 * (1 - 6*phi + 6 * phi^2)
    int dof = el.GetDof();
    int dim = el.GetDim();
    shape.SetSize(dof);
    dshape.SetSize(dof, dim);
    invdfdx.SetSize(dim);
    vec.SetSize(dim);
    pointflux.SetSize(dim);

    elvect.SetSize(dof);
    elvect = 0.0;

    const mfem::IntegrationRule* ir =
        &mfem::IntRules.Get(el.GetGeomType(), 2 * el.GetOrder() + Tr.OrderW());

    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      el.CalcShape(ip, shape);
      Tr.SetIntPoint(&ip);

      // Given u, compute (w'(u)+calphad, v), v is shape function

      auto mobi = (*mob).Eval(Tr, ip);
      auto lamb = (*lambda).Eval(Tr, ip);
      auto _omega = (*omega).Eval(Tr, ip);
      auto phi = (elfun * shape);
      auto wp = _omega * firstDerivativedoubleWellPotential(phi);
      auto calphadContribution = (*Calphad).Eval(Tr, ip);

      double fun_val = mobi * (wp + calphadContribution);
      double w = ip.weight * Tr.Weight() * fun_val;
      add(elvect, w, shape, elvect);

      // Given u, compute (grad(u), grad(v)), v is shape function. Ref:
      // DiffusionIntegrator::AssembleElementVector()
      CalcAdjugate(Tr.Jacobian(), invdfdx);
      dshape.MultTranspose(elfun, vec);
      invdfdx.MultTranspose(vec, pointflux);
      double ww = ip.weight / Tr.Weight();
      pointflux *= ww;
      pointflux *= mobi * lamb;
      invdfdx.Mult(pointflux, vec);
      dshape.AddMult(vec, elvect);
    }
  };
  virtual void AssembleElementGrad(const mfem::FiniteElement& el,
                                   mfem::ElementTransformation& Tr,
                                   const mfem::Vector& elfun,
                                   mfem::DenseMatrix& elmat) {
    // Linearized DoubleWell contribution
    // W' = W'(phi_n)+W''(phi_n)*(phi_n+1 - phi_n)
    // with :
    // W = phi^2 * (1-phi)^2
    // W' = 2 phi * (1-phi) * (1-2*phi)
    // w'' = 2 * (1 - 6*phi + 6 * phi^2)
    int dof = el.GetDof();
    int dim = el.GetDim();
    dshapedxt.SetSize(dof, dim);
    dshape.SetSize(dof, dim);
    shape.SetSize(dof);
    elmat.SetSize(dof);
    elmat = 0.0;

    const mfem::IntegrationRule* ir =
        &mfem::IntRules.Get(el.GetGeomType(), 2 * el.GetOrder() + Tr.OrderW());

    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      el.CalcShape(ip, shape);
      el.CalcDShape(ip, dshape);
      Tr.SetIntPoint(&ip);

      auto mobi = (*mob).Eval(Tr, ip);
      auto lamb = (*lambda).Eval(Tr, ip);
      auto _omega = (*omega).Eval(Tr, ip);
      auto phi = (elfun * shape);
      auto wpp = _omega * secondDerivativedoubleWellPotential(phi);

      // Compute (grad(du), grad(v)).  Ref:
      // DiffusionIntegrator::AssembleElementMatrix()
      double w = ip.weight / Tr.Weight();
      w *= _omega * lamb;
      Mult(dshape, Tr.AdjugateJacobian(), dshapedxt);  //
      AddMult_a_AAt(w, dshapedxt, elmat);

      // TODO rajouter le terme en hp(phi)
      // Compute w''(u)*(du,v), v is shape function
      double fun_val = _omega * wpp * ip.weight * Tr.Weight();  // w''(u)
      AddMult_a_VVt(fun_val, shape, elmat);  // w''(u)*(du, v)
    }
  };

  virtual ~AllenCahnSpecializedNLFormIntegrator(){};
};
#endif /* SOLVERS_ALLENCAHNSPECIALIZEDNLFORMINTEGRATOR_H_ \
          */

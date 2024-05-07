/*
 * Class for integrating the Nonlinear form a(u,v) := (Q . Qnl(u) . grad(u), grad(v)) in
 *    either 2D or 3D and where Q is an optional coefficient (of type scalar, matrix or
 *    diagonal matrix), Qnl(u) is an nonlinear coefficient w.r.t. u and u and v are in
 *    H1.
 */
#include "mfem.hpp"

#ifndef SOLVERS_DIFFUSIONNLFINTEGRATOR_HPP_
#define SOLVERS_DIFFUSIONNLFINTEGRATOR_HPP_
/** Function representing a nonlinear coefficient (density) for the given state
 *  (pressure). Used in DiffusionNLFIntegrator::AssembleElementGrad. */
class NonlinearCoefficient : public mfem::GridFunctionCoefficient {
 private:
  // variables members
  double rho0, beta, u0;

 public:
  NonlinearCoefficient(double rho_, double bt_, double p0)
      : mfem::GridFunctionCoefficient(), rho0(rho_), beta(bt_), u0(p0) {}

  NonlinearCoefficient(mfem::GridFunction *u_, double rho_, double bt_, double p0)
      : mfem::GridFunctionCoefficient(u_), rho0(rho_), beta(bt_), u0(p0) {}

  virtual double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) {
    double u = mfem::GridFunctionCoefficient::Eval(T, ip);
    return rho0 * std::exp(beta * (u - u0));
  }

  virtual double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip,
                      const double &u) {
    return rho0 * std::exp(beta * (u - u0));
  }

  virtual double EvalGrad(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip,
                          const double &u) {
    // Set the size of the vector to match the problem dimensionality
    // double rho = Eval(T, ip);
    return rho0 * beta * std::exp(beta * (u - u0));
    ;
  }

  virtual void Read(std::istream &in) {}

  virtual ~NonlinearCoefficient() {}
};

class DiffusionNLFIntegrator : public mfem::NonlinearFormIntegrator {
 protected:
  mfem::ConstantCoefficient *Q;
  mfem::VectorCoefficient *VQ;
  mfem::MatrixCoefficient *MQ;
  mfem::SymmetricMatrixCoefficient *SMQ;
  NonlinearCoefficient *QNL;

 private:
  mfem::Vector vec, vecdxt, pointflux, shape;
#ifndef MFEM_THREAD_SAFE
  mfem::DenseMatrix dshape, dshapedxt, Jadj, M, dshapedxt_m;
  mfem::Vector D, buf;
#endif

 public:
  /// Construct a diffusion nonliner integrator with coefficient Q = 1
  DiffusionNLFIntegrator() : Q(NULL), VQ(NULL), MQ(NULL), SMQ(NULL), QNL(NULL) {}

  /// Construct a diffusion integrator with a scalar coefficient q and nonlinear
  ///   coefficient qq
  DiffusionNLFIntegrator(mfem::ConstantCoefficient &q)
      : Q(&q), VQ(NULL), MQ(NULL), SMQ(NULL), QNL(NULL) {}

  /// Construct a diffusion integrator with a scalar coefficient q and nonlinear
  ///   coefficient qq
  DiffusionNLFIntegrator(mfem::ConstantCoefficient &q, NonlinearCoefficient &qq)
      : Q(&q), VQ(NULL), MQ(NULL), SMQ(NULL), QNL(&qq) {}

  /// Construct a diffusion integrator with a vector coefficient q and nonlinear
  ///   coefficient qq
  DiffusionNLFIntegrator(mfem::VectorCoefficient &q, NonlinearCoefficient &qq)
      : Q(NULL), VQ(&q), MQ(NULL), SMQ(NULL), QNL(&qq) {}

  /// Construct a diffusion integrator with a matrix coefficient q and nonlinear
  ///   coefficient qq
  DiffusionNLFIntegrator(mfem::MatrixCoefficient &q, NonlinearCoefficient &qq)
      : Q(NULL), VQ(NULL), MQ(&q), SMQ(NULL), QNL(&qq) {}

  /// Construct a diffusion integrator with a symmetric matrix coefficient q and
  ///   nonlinear coefficient qq
  DiffusionNLFIntegrator(mfem::SymmetricMatrixCoefficient &q, NonlinearCoefficient &qq)
      : Q(NULL), VQ(NULL), MQ(NULL), SMQ(&q), QNL(&qq) {}

  static const mfem::IntegrationRule &GetRule(const mfem::FiniteElement &fe,
                                              mfem::ElementTransformation &T);

  /// Given a elfun values perform the local action of the NonlinearFormIntegrator
  virtual void AssembleElementVector(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr,
                                     const mfem::Vector &elfun, mfem::Vector &elvect);

  /// Assemble the local gradient matrix
  virtual void AssembleElementGrad(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr,
                                   const mfem::Vector &elfun, mfem::DenseMatrix &elmat);
};

////////////////////////////////////////////////////////////////////////////////
//////////////////////// DiffusionNLFIntegrator methods ////////////////////////
////////////////////////////////////////////////////////////////////////////////

// implementation based on:
//  - Ref. mfem/fem/nonlininteg.cpp::VectorConvectionNLFIntegrator::GetRule
//
const mfem::IntegrationRule &DiffusionNLFIntegrator::GetRule(const mfem::FiniteElement &fe,
                                                             mfem::ElementTransformation &T) {
  const int order = 2 * fe.GetOrder() + T.OrderGrad(&fe);
  return mfem::IntRules.Get(fe.GetGeomType(), order);
}

// implementation based on:
//  - Ref. mfem/fem/bilininteg.cpp::DiffusionIntegrator::AssembleElementVector
//  - Ref. mfem/fem/nonlininteg.cpp::VectorConvectionNLFIntegrator::AssembleElementVector
//
// compute (Q * Qnl(u) * grad(u), grad(v))
void DiffusionNLFIntegrator::AssembleElementVector(const mfem::FiniteElement &el,
                                                   mfem::ElementTransformation &Tr,
                                                   const mfem::Vector &elfun,
                                                   mfem::Vector &elvect) {
  // local el. info.
  const int dof_per_el = el.GetDof();
  const int dim = el.GetDim();
  int spaceDim = Tr.GetSpaceDim();

  // verify coefficient type in DiffusionNLFIntegrator
  if (VQ) {
    MFEM_VERIFY(VQ->GetVDim() == spaceDim, "Unexpected dimension for VectorCoefficient");
  }
  if (MQ) {
    MFEM_VERIFY(MQ->GetWidth() == spaceDim, "Unexpected width for MatrixCoefficient");
    MFEM_VERIFY(MQ->GetHeight() == spaceDim, "Unexpected height for MatrixCoefficient");
  }

#ifdef MFEM_THREAD_SAFE
  DenseMatrix dshape(dof_per_el, dim), M(MQ ? spaceDim : 0), Jadj(dim, spaceDim);
  Vector D(VQ ? VQ->GetVDim() : 0);
#else
  shape.SetSize(dof_per_el);
  dshape.SetSize(dof_per_el, dim);
  Jadj.SetSize(dim, spaceDim);
  M.SetSize(MQ ? spaceDim : 0);
  D.SetSize(VQ ? VQ->GetVDim() : 0);
#endif
  vec.SetSize(dim);
  vecdxt.SetSize((VQ || MQ) ? spaceDim : 0);
  pointflux.SetSize(spaceDim);

  elvect.SetSize(dof_per_el);
  elvect = 0.0;
  double w, rho;

  // integration rule
  const mfem::IntegrationRule *ir = IntRule ? IntRule : &GetRule(el, Tr);

  // loop over quad./integration point
  for (int i = 0; i < ir->GetNPoints(); i++) {
    const mfem::IntegrationPoint &ip = ir->IntPoint(i);
    el.CalcShape(ip, shape);    // phi
    el.CalcDShape(ip, dshape);  // dphi

    Tr.SetIntPoint(&ip);                // fixes ip (used in the next function)
    CalcAdjugate(Tr.Jacobian(), Jadj);  // Jadj = adj(J(xq)) = detJ*J^{-1}
    w = ip.weight / Tr.Weight();        // wq/det(J)

    // density(u)
    double u = shape * elfun;
    rho = 1.0;
    if (QNL) {
      rho = QNL->Eval(Tr, ip, u);  //=rho(xq)
    }

    // Tensor/Stiff coefficient
    if (!MQ && !VQ) {
      dshape.MultTranspose(elfun, vec);    //= dphi^t*elfun = dphi^t*u
      Jadj.MultTranspose(vec, pointflux);  //= Jadj^t*dphi^t*u
      if (Q) {
        w *= Q->Eval(Tr, ip);  // wq/det(J)*Q(xq)
      }
    } else {
      dshape.MultTranspose(elfun, vec);  // vec = dphi^t*elfun = dphi^t*u
      Jadj.MultTranspose(vec, vecdxt);   // Jadj^t*dphi^t*u
      if (MQ) {
        MQ->Eval(M, Tr, ip);        // Q(xq)
        M.Mult(vecdxt, pointflux);  //= Q(xq)*Jadj^t*dphi^t*u
      } else {
        VQ->Eval(D, Tr, ip);  // Q(xq)
        for (int j = 0; j < spaceDim; ++j) {
          pointflux[j] = D[j] * vecdxt[j];  // Q[i]*grad(u)[i]
        }
      }
    }
    pointflux *= (w * rho);       //= wq/det(J)*rho(xq)*Q(xq)*Jadj^t*dphi*u
    Jadj.Mult(pointflux, vec);    //= wq*Jadj^t*rho(xq)*Q(xq)*Jadj^t*dphi^t*u
    dshape.AddMult(vec, elvect);  // wq*[Jadj^t*dphi_i]*rho(xq)*Q(xq)*[Jadj^t*dphi_j^t]*u
  }
}

// implementation based on:
//  - Ref. mfem/fem/bilininteg.cpp::DiffusionIntegrator::AssembleElementVector
//  - Ref. mfem/fem/bilininteg.cpp::MassIntegrator::AssembleElementMatrix2
//  - Ref. mfem/fem/nonlininteg.cpp::VectorConvectionNLFIntegrator::AssembleElementVector
//
// compute jac. of (Q * Qnl(u) * grad(u), grad(v))
//      (Q * Qnl(u) * grad(du), grad(v)) + (Qnl'(u) * Q * grad(u) * du,grad(v))
void DiffusionNLFIntegrator::AssembleElementGrad(const mfem::FiniteElement &el,
                                                 mfem::ElementTransformation &Tr,
                                                 const mfem::Vector &elfun,
                                                 mfem::DenseMatrix &elmat) {
  // local el. info.
  int dof_per_el = el.GetDof();
  int dim = el.GetDim();
  int spaceDim = Tr.GetSpaceDim();
  bool square = (dim == spaceDim);

  // verify coefficient type in DiffusionNLFIntegrator
  if (VQ) {
    MFEM_VERIFY(VQ->GetVDim() == spaceDim, "Unexpected dimension for VectorCoefficient");
  }
  if (MQ) {
    MFEM_VERIFY(MQ->GetWidth() == spaceDim, "Unexpected width for MatrixCoefficient");
    MFEM_VERIFY(MQ->GetHeight() == spaceDim, "Unexpected height for MatrixCoefficient");
  }

#ifdef MFEM_THREAD_SAFE
  mfem::DenseMatrix dshape(dof_per_el, dim), dshapedxt(dof_per_el, spaceDim);
  mfem::DenseMatrix dshapedxt_m(dof_per_el, MQ ? spaceDim : 0), Jadj(dim, spaceDim);
  mfem::Vector D(VQ ? VQ->GetVDim() : 0);
  mfem::Vector buf(dof_per_el), vec(dim);
#else
  shape.SetSize(dof_per_el);        // phi
  dshape.SetSize(dof_per_el, dim);  // dphi
  Jadj.SetSize(dim, spaceDim);      // det(J)*J^{-1}
  dshapedxt.SetSize(dof_per_el, spaceDim);
  dshapedxt_m.SetSize(dof_per_el, MQ ? spaceDim : 0);
  M.SetSize(MQ ? spaceDim : 0);
  D.SetSize(VQ ? VQ->GetVDim() : 0);
  buf.SetSize(dof_per_el);
#endif
  vec.SetSize(dim);
  vecdxt.SetSize((VQ || MQ) ? spaceDim : 0);
  pointflux.SetSize(spaceDim);
  elmat.SetSize(dof_per_el);
  elmat = 0.0;
  double w, rho, drho;

  // integration rule
  const mfem::IntegrationRule *ir = IntRule ? IntRule : &GetRule(el, Tr);

  // loop over quad./integration point
  for (int i = 0; i < ir->GetNPoints(); i++) {
    const mfem::IntegrationPoint &ip = ir->IntPoint(i);
    el.CalcShape(ip, shape);    // phi = [N1 N2 ... Nn]
                                //        |dN1dx dN1dy dN1dz|
    el.CalcDShape(ip, dshape);  // dphi = |dN2dx dN2dy dN2dz|
                                //        |dNidx dNidy dNidz|
                                //        | ...   ...   ... |
                                //        |dNndx dNndy dNndz|

    Tr.SetIntPoint(&ip);  // fixes ip (used in the next function)
    w = Tr.Weight();      // wq
    w = ip.weight / (square ? w : w * w * w);
    // AdjugateJacobian = / adj(J),         if J is square
    //                    \ adj(J^t.J).J^t, otherwise
    // adj(J(xq)) = detJ*J^{-1}
    CalcAdjugate(Tr.Jacobian(), Jadj);  // Jadj = adj(J(xq)) = detJ*J^{-1}
    Mult(dshape, Jadj, dshapedxt);      //=detJ*J^{-1}*dphi

    rho = 1.0;
    drho = 0.0;
    double u = shape * elfun;
    if (QNL) {
      rho = QNL->Eval(Tr, ip, u);       //=rho(xq)
      drho = QNL->EvalGrad(Tr, ip, u);  //=drho(xq)
    }

    // ====================================
    // compute (rho*kappa*grad(du),grad(v))
    if (MQ) {
      MQ->Eval(M, Tr, ip);              // Q(xq)
      M *= w * rho;                     //=wq/det(J)*Q(xq)*rho(xq)
      Mult(dshapedxt, M, dshapedxt_m);  //=wq/det(J)*Q(xq)*detJ*J^{-1}*dphi
      AddMultABt(dshapedxt_m, dshapedxt,
                 elmat);  //=wq/det(J)*Q(xq)*detJ*J^{-1}*dphi*[detJ*J^{-1}*dphi]^t
    } else if (VQ) {
      VQ->Eval(D, Tr, ip);               // Q(xq)
      D *= w * rho;                      //=wq/det(J)*Q(xq)*rho(xq)
      AddMultADAt(dshapedxt, D, elmat);  //=wq/det(J)*Q(xq)*detJ*J^{-1}*dphi*[detJ*J^{-1}*dphi]^t
    } else {
      if (Q) {
        w *= Q->Eval(Tr, ip) * rho;  //=wq/det(J)*Q(xq)*rho(xq)
      }
      AddMult_a_AAt(w, dshapedxt, elmat);  //=wq/det(J)*Q(xq)*detJ*J^{-1}*dphi*[detJ*J^{-1}*dphi]^t
    }

    // ====================================
    // compute drho/du*kappa*grad(u)*(du,grad(v))
    dshape.MultTranspose(elfun, vec);  //= dphi^t*elfun = dphi^t*u
    if (MQ) {
      Jadj.MultTranspose(vec, vecdxt);  //= Jadj^t*dphi^t*u
      MQ->Eval(M, Tr, ip);              // Q(xq)
      M.Mult(vecdxt, pointflux);        //=Q(xq)*Jadj^t*dphi^t*u
      dshapedxt.Mult(pointflux, buf);   //=Jadj*dphi*Q(xq)*Jadj^t*dphi^t*u
      buf *= (w * drho);                //=wq*invJ*dphi*Q(xq)*drho*invJ^t*dphi^t*u*det(J)
      AddMultVWt(buf, shape, elmat);
    } else if (VQ) {
      Jadj.MultTranspose(vec, vecdxt);  //= Jadj^t*dphi^t*u
      VQ->Eval(D, Tr, ip);              // Q(xq)
      for (int j = 0; j < spaceDim; ++j) {
        //=wq*invJ*dphi*Q(xq)*drho*invJ^t*dphi^t*u*det(J)
        pointflux[j] = w * drho * D[j] * vecdxt[j];
      }
      dshapedxt.Mult(pointflux, buf);  //=Jadj*dphi*Q(xq)*drho*Jadj^t*dphi^t*u
      AddMultVWt(buf, shape, elmat);
    } else {
      if (Q) {
        w *= Q->Eval(Tr, ip) * drho;  //=wq/det(J)*Q(xq)
      }
      pointflux *= w;                  //=wq*invJ*dphi*Q(xq)*drho*invJ^t*dphi^t*u*det(J)
      dshapedxt.Mult(pointflux, buf);  //=Jadj*dphi*Q(xq)*Jadj^t*dphi^t*u
      AddMultVWt(buf, shape, elmat);
    }
  }
}

#endif /* SOLVERS_DIFFUSIONNLFINTEGRATOR_HPP_ */
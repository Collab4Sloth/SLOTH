/*
 * Copyright © CEA 2022
 *
 * \brief Conduction operator used for test and learn
 *
 * \file ConductionOpertor.hpp
 * \author ci230846 (from MFEM example)
 * \date 02/12/2022
 */

#include <map>
#include <memory>
#include "mfem.hpp"

#pragma once
/*
 *  Class ConductionOperator represents the right-hand side of the above ODE.
 */
class ConductionOperator : public mfem::TimeDependentOperator {
 protected:
  mfem::FiniteElementSpace &fespace;
  mfem::Array<int> ess_tdof_list;  // this list remains empty for pure Neumann b.c.

  mfem::BilinearForm *M;
  mfem::BilinearForm *K;

  mfem::SparseMatrix Mmat, Kmat;
  mfem::SparseMatrix *T;  // T = M + dt K
  double current_dt;

  mfem::CGSolver M_solver;  // Krylov solver for inverting the mass matrix M
  mfem::DSmoother M_prec;   // Preconditioner for the mass matrix M

  // mfem::CGSolver T_solver;  // Implicit solver for T = M + dt K
  // mfem::DSmoother T_prec;   // Preconditioner for the implicit solver
  mfem::UMFPackSolver T_solver;  // Implicit solver for T = M + dt K
  // mfem::Solver T_prec;    // Preconditioner for the implicit solver

  double alpha, kappa;

  mutable mfem::Vector z;  // auxiliary vector

 public:
  ConductionOperator(mfem::FiniteElementSpace &f, double alpha, double kappa, mfem::Vector &u);

  virtual void Mult(const mfem::Vector &u, mfem::Vector &du_dt) const;
  /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
      This is the only requirement for high-order SDIRK implicit integration.*/
  virtual void ImplicitSolve(const double dt, const mfem::Vector &u, mfem::Vector &k);

  /// Update the diffusion BilinearForm K using the given true-dof vector `u`.
  void SetParameters(const mfem::Vector &u);

  virtual ~ConductionOperator();
};

ConductionOperator::ConductionOperator(mfem::FiniteElementSpace &f, double al, double kap,
                                       mfem::Vector &u)
    : mfem::TimeDependentOperator(f.GetTrueVSize(), 0.0),
      fespace(f),
      ess_tdof_list(),
      M(NULL),
      K(NULL),
      T(NULL),
      current_dt(0.0),
      z(height) {
  const double rel_tol = 1e-8;
  const double abs_tol = 1e-8;

  ///////////////////////////
  // BC
  ///////////////////////////
  mfem::Array<int> Dirichlet_bdr_1(f.GetMesh()->bdr_attributes.Max());
  mfem::Array<int> Dirichlet_bdr_3(f.GetMesh()->bdr_attributes.Max());
  mfem::Array<int> Dirichlet_bdr(f.GetMesh()->bdr_attributes.Max());
  Dirichlet_bdr_1 = 0;
  Dirichlet_bdr_3 = 0;
  Dirichlet_bdr[0] = 0;
  Dirichlet_bdr[1] = 0;
  Dirichlet_bdr[2] = 0;
  Dirichlet_bdr[3] = 0;
  Dirichlet_bdr_1[1] = 0;
  Dirichlet_bdr_3[3] = 0;

  // fespace.GetEssentialTrueDofs(Dirichlet_bdr_1, ess_tdof_list);
  // static constexpr double ub1 = 1.0;
  // u.SetSubVector(ess_tdof_list, ub1);
  // fespace.GetEssentialTrueDofs(Dirichlet_bdr_3, ess_tdof_list);
  // static constexpr double ub3 = 3.0;
  // u.SetSubVector(ess_tdof_list, ub3);
  fespace.GetEssentialTrueDofs(Dirichlet_bdr, ess_tdof_list);

  M = new mfem::BilinearForm(&fespace);
  M->AddDomainIntegrator(new mfem::MassIntegrator());
  M->Assemble();
  M->FormSystemMatrix(ess_tdof_list, Mmat);

  M_solver.iterative_mode = false;
  M_solver.SetRelTol(rel_tol);
  M_solver.SetAbsTol(abs_tol);
  M_solver.SetMaxIter(30);
  M_solver.SetPrintLevel(0);
  M_solver.SetPreconditioner(M_prec);
  M_solver.SetOperator(Mmat);

  alpha = al;
  kappa = kap;

  // T_solver.iterative_mode = false;
  // T_solver.SetRelTol(rel_tol);
  // T_solver.SetAbsTol(abs_tol);
  // T_solver.SetMaxIter(100);
  T_solver.SetPrintLevel(1);
  // T_solver.SetPreconditioner(T_prec);

  // T_solver = new mfem::UMFPackSolver;
  // T_prec = null;

  SetParameters(u);
}

void ConductionOperator::Mult(const mfem::Vector &u, mfem::Vector &du_dt) const {
  // Compute:
  //    du_dt = M^{-1}*-Ku
  // for du_dt, where K is linearized by using u from the previous timestep
  // z = K.u
  Kmat.Mult(u, z);
  z.Neg();  // z = -z
  // du_dt = M * K.u

  // for (int i = 0; i < ess_tdof_list.Size(); i++)
  //   z(ess_tdof_list[i]) = 0.0;  // set Dirichlet condition by hand

  M_solver.Mult(z, du_dt);
}

void ConductionOperator::ImplicitSolve(const double dt, const mfem::Vector &u,
                                       mfem::Vector &du_dt) {
  // Solve the equation:
  //    du_dt = M^{-1}*[-K(u + dt*du_dt)]
  // for du_dt, where K is linearized by using u from the previous timestep
  if (!T) {
    T = Add(1.0, Mmat, dt, Kmat);
    current_dt = dt;
    T_solver.SetOperator(*T);
  }
  MFEM_VERIFY(dt == current_dt, "");  // SDIRK methods use the same dt
  Kmat.Mult(u, z);
  z.Neg();

  T_solver.Mult(z, du_dt);
  du_dt.SetSubVector(ess_tdof_list, 0.0);
}

void ConductionOperator::SetParameters(const mfem::Vector &u) {
  mfem::GridFunction u_alpha_gf(&fespace);
  u_alpha_gf.SetFromTrueDofs(u);
  for (int i = 0; i < u_alpha_gf.Size(); i++) {
    u_alpha_gf(i) = kappa + alpha * u_alpha_gf(i);
  }

  delete K;
  K = new mfem::BilinearForm(&fespace);

  mfem::GridFunctionCoefficient u_coeff(&u_alpha_gf);

  K->AddDomainIntegrator(new mfem::DiffusionIntegrator(u_coeff));
  K->Assemble();
  K->FormSystemMatrix(ess_tdof_list, Kmat);
  delete T;
  T = NULL;  // re-compute T on the next ImplicitSolve
}

ConductionOperator::~ConductionOperator() {
  delete T;
  delete M;
  delete K;
}

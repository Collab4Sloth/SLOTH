/*
 * Copyright © CEA 2022
 *
 * \brief PhaseField time-dependent operator built on the basis of ex16.cpp from MFEM repository
 *
 *   After spatial discretization, the phasefield model can be written as:
 *      dphi/dt = M^{-1}(-Kphi)
 *   where phi denotes the phase-field variable and K the phase-field operator (that does not
 *   depends on time)
 *
 * \file ReducedOperator.hpp
 * \author ci230846
 * \date 20/01/2022
 */

#include "mfem.hpp"

#pragma once
/*
 *  Class PhaseFieldReducedOperator
 */
class PhaseFieldReducedOperator : public mfem::Operator {
 private:
  // Mass matrix
  mfem::BilinearForm *M;
  // PhaseField Matrix
  mfem::NonlinearForm *N;
  // Jacobian matrix
  mutable mfem::SparseMatrix *Jacobian;

  // Time step
  double dt;
  // Unknown
  const mfem::Vector *unk;
  mutable mfem::Vector z;

 public:
  PhaseFieldReducedOperator(mfem::BilinearForm *M_, mfem::NonlinearForm *N_);

  /// Set current dt, unk values - needed to compute action and Jacobian.
  void SetParameters(double dt_, const mfem::Vector *unk_);

  /// Compute y = N(unk + dt*k) + M k
  void Mult(const mfem::Vector &k, mfem::Vector &y) const;

  /// Compute y = dt*grad_N(unk + dt*k) + M
  mfem::Operator &GetGradient(const mfem::Vector &k) const;
  ~PhaseFieldReducedOperator();
};

PhaseFieldReducedOperator::PhaseFieldReducedOperator(mfem::BilinearForm *M_,
                                                     mfem::NonlinearForm *N_)
    : Operator(N_->Height()), M(M_), N(N_), Jacobian(NULL), dt(0.0), unk(NULL), z(height) {}

/// Set current dt, unk values - needed to compute action and Jacobian.
void PhaseFieldReducedOperator::SetParameters(double dt_, const mfem::Vector *unk_) {
  dt = dt_;
  unk = unk_;
}

/// Compute y = N(unk + dt*k) + M k
void PhaseFieldReducedOperator::Mult(const mfem::Vector &k, mfem::Vector &y) const {
  add(*unk, dt, k, z);
  N->Mult(z, y);
  M->AddMult(k, y);
}

mfem::Operator &PhaseFieldReducedOperator::GetGradient(const mfem::Vector &k) const {
  // std::cout << " PhaseFieldReducedOperator::GetGradient  " << std::endl;
  delete Jacobian;
  Jacobian = Add(1.0, M->SpMat(), 0.0, M->SpMat());
  add(*unk, dt, k, z);
  mfem::SparseMatrix *grad_N = dynamic_cast<mfem::SparseMatrix *>(&N->GetGradient(z));
  Jacobian->Add(dt, *grad_N);
  return *Jacobian;
}

PhaseFieldReducedOperator::~PhaseFieldReducedOperator() { delete Jacobian; }

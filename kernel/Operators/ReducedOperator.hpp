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

#include "mfem.hpp" // NOLINT [no include the directory when naming mfem include file]

#pragma once
/*
 *  Class PhaseFieldReducedOperator
 */
class PhaseFieldReducedOperator : public mfem::Operator {
 private:
  // Mass matrix
  mfem::BilinearForm *M_;
  // PhaseField Matrix
  mfem::NonlinearForm *N_;
  // Jacobian matrix
  mutable mfem::SparseMatrix *Jacobian;

  // Time step
  double dt_;
  // Unknown
  const mfem::Vector *unk_;
  mutable mfem::Vector z;

 public:
  PhaseFieldReducedOperator(mfem::BilinearForm *M, mfem::NonlinearForm *N);

  /// Set current dt, unk values - needed to compute action and Jacobian.
  void SetParameters(double dt, const mfem::Vector *unk);

  /// Compute y = N(unk + dt*k) + M k
  void Mult(const mfem::Vector &k, mfem::Vector &y) const;

  /// Compute y = dt*grad_N(unk + dt*k) + M
  mfem::Operator &GetGradient(const mfem::Vector &k) const;
  ~PhaseFieldReducedOperator();
};

PhaseFieldReducedOperator::PhaseFieldReducedOperator(mfem::BilinearForm *M, mfem::NonlinearForm *N)
    : Operator(N->Height()), M_(M), N_(N), Jacobian(NULL), dt_(0.0), unk_(NULL), z(height) {}

/// Set current dt, unk values - needed to compute action and Jacobian.
void PhaseFieldReducedOperator::SetParameters(double dt, const mfem::Vector *unk) {
  dt_ = dt;
  unk_ = unk;
}

/// Compute y = N(unk + dt*k) + M k
void PhaseFieldReducedOperator::Mult(const mfem::Vector &k, mfem::Vector &y) const {
  add(*unk_, dt_, k, z);
  N_->Mult(z, y);
  M_->AddMult(k, y);
}

mfem::Operator &PhaseFieldReducedOperator::GetGradient(const mfem::Vector &k) const {
  delete Jacobian;
  Jacobian = Add(1.0, M_->SpMat(), 0.0, M_->SpMat());
  add(*unk_, dt_, k, z);
  mfem::SparseMatrix *grad_N = dynamic_cast<mfem::SparseMatrix *>(&N_->GetGradient(z));
  Jacobian->Add(dt_, *grad_N);
  return *Jacobian;
}

PhaseFieldReducedOperator::~PhaseFieldReducedOperator() { delete Jacobian; }

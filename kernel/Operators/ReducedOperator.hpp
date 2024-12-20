/**
 * @file ReducedOperator.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief
 * @version 0.1
 * @date 2024-07-30
 *
 * Copyright CEA (c) 2024
 *
 */
#include <memory>

#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
#pragma once
/*
 *  Class PhaseFieldReducedOperator
 */
class PhaseFieldReducedOperator : public mfem::Operator {
 private:
  // Mass matrix
  mfem::ParBilinearForm *M_;
  // PhaseField Matrix
  mfem::ParNonlinearForm *N_;
  // Jacobian matrix
  mutable mfem::HypreParMatrix *Jacobian;

  // Time step
  double dt_;
  // Unknown
  const mfem::Vector *unk_;
  mutable mfem::Vector z;

  const mfem::Array<int> &ess_tdof_list;

 public:
  PhaseFieldReducedOperator(mfem::ParBilinearForm *M, mfem::ParNonlinearForm *N,
                            const mfem::Array<int> &ess_tdof);

  /// Set current dt, unk values - needed to compute action and Jacobian.
  void SetParameters(double dt, const mfem::Vector *unk);

  /// Compute y = N(unk + dt*k) + M k
  void Mult(const mfem::Vector &k, mfem::Vector &y) const;

  /// Compute y = dt*grad_N(unk + dt*k) + M
  mfem::Operator &GetGradient(const mfem::Vector &k) const;
  ~PhaseFieldReducedOperator();
};

/**
 * @brief Construct a new Phase Field Reduced Operator:: Phase Field Reduced Operator object
 *
 * @param M
 * @param N
 */
PhaseFieldReducedOperator::PhaseFieldReducedOperator(mfem::ParBilinearForm *M,
                                                     mfem::ParNonlinearForm *N,
                                                     const mfem::Array<int> &ess_tdof)
    : Operator(N->ParFESpace()->TrueVSize()),
      M_(M),
      N_(N),
      Jacobian(NULL),
      dt_(0.0),
      unk_(NULL),
      z(height),
      ess_tdof_list(ess_tdof) {}

/**
 * @brief  Set current dt, unk values - needed to compute action and Jacobian.
 *
 * @param dt
 * @param unk
 */
void PhaseFieldReducedOperator::SetParameters(double dt, const mfem::Vector *unk) {
  dt_ = dt;
  unk_ = unk;
}

/**
 * @brief  Compute y = N(unk + dt*k) + M k
 *
 * @param k
 * @param y
 */
void PhaseFieldReducedOperator::Mult(const mfem::Vector &k, mfem::Vector &y) const {
  add(*unk_, dt_, k, z);
  N_->Mult(z, y);
  M_->TrueAddMult(k, y);
  y.SetSubVector(ess_tdof_list, 0.0);
}

/**
 * @brief  Compute Jacobian
 *
 * @param k
 * @return mfem::Operator&
 */
mfem::Operator &PhaseFieldReducedOperator::GetGradient(const mfem::Vector &k) const {
  if (Jacobian != nullptr) {
    delete Jacobian;
  }

  std::unique_ptr<mfem::SparseMatrix> localJ(Add(1.0, M_->SpMat(), 0.0, M_->SpMat()));

  add(*unk_, dt_, k, z);

  localJ->Add(dt_, N_->GetLocalGradient(z));
  Jacobian = M_->ParallelAssemble(localJ.get());

  std::unique_ptr<mfem::HypreParMatrix> Je(Jacobian->EliminateRowsCols(ess_tdof_list));

  return *Jacobian;
}

/**
 * @brief Destroy the Phase Field Reduced Operator:: Phase Field Reduced Operator object
 *
 */
PhaseFieldReducedOperator::~PhaseFieldReducedOperator() { delete Jacobian; }

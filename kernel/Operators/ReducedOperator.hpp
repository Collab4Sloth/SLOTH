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
  mfem::ParBlockNonlinearForm *LHS_;
  // PhaseField Matrix
  mfem::ParBlockNonlinearForm *N_;
  // Jacobian matrix
  mutable mfem::HypreParMatrix *Jacobian;

  // Time step
  double dt_;
  // Unknown
  const mfem::Vector *unk_;
  mutable mfem::Vector z;

  const std::vector<mfem::Array<int>> &ess_tdof_list;

 public:
  PhaseFieldReducedOperator(mfem::ParBlockNonlinearForm *M, mfem::ParBlockNonlinearForm *N,
                            const std::vector<mfem::Array<int>> &ess_tdof);

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
PhaseFieldReducedOperator::PhaseFieldReducedOperator(mfem::ParBlockNonlinearForm *LHS,
                                                     mfem::ParBlockNonlinearForm *N,
                                                     const std::vector<mfem::Array<int>> &ess_tdof)
    : Operator(N->Height()),
      // : Operator(N->ParFESpace()->TrueVSize()),
      LHS_(LHS),
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
  // M_->TrueAddMult(k, y);
  LHS_->AddMult(k, y);
  // TODO(cci) adapt
  // y.SetSubVector(ess_tdof_list, 0.0);
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
  add(*unk_, dt_, k, z);
  const mfem::Array<int> offsets = N_->GetBlockOffsets();
  const int fes_size = offsets.Size() - 1;
  // Gets gradients of N_ and LHS_
  mfem::Operator &LHS_grad = LHS_->GetGradient(z);
  mfem::Operator &N_grad = N_->GetGradient(z);
  // Converts operators into BlockOperator
  mfem::BlockOperator *LHS_block_grad = dynamic_cast<mfem::BlockOperator *>(&LHS_grad);
  mfem::BlockOperator *N_block_grad = dynamic_cast<mfem::BlockOperator *>(&N_grad);
  mfem::Array2D<mfem::HypreParMatrix *> tmp_blocks(fes_size, fes_size);
  std::vector<mfem::HypreParMatrix *> blocks_to_delete;

  for (int i = 0; i < fes_size; ++i) {
    for (int j = 0; j < fes_size; ++j) {
      mfem::Operator *LHS_block = &(LHS_block_grad->GetBlock(i, j));
      mfem::Operator *N_block = &(N_block_grad->GetBlock(i, j));

      mfem::HypreParMatrix *LHS_sparse_block = dynamic_cast<mfem::HypreParMatrix *>(LHS_block);
      mfem::HypreParMatrix *N_sparse_block = dynamic_cast<mfem::HypreParMatrix *>(N_block);

      if (LHS_sparse_block && N_sparse_block) {
        // mfem::HypreParMatrix *block = new mfem::HypreParMatrix(*LHS_sparse_block);
        // std::cout << " coucou2" << std::endl;
        // block->Add(dt_, *N_sparse_block);
        // std::cout << " coucou3" << std::endl;
        // // Jacobian->SetBlock(i, j, block);
        // tmp_blocks(i, j) = block;

        //
        // TODO for ch
        if (i == 1) {
          mfem::HypreParMatrix *block = new mfem::HypreParMatrix(*LHS_sparse_block);
          block->Add(dt_, *N_sparse_block);
          tmp_blocks(i, j) = block;
          blocks_to_delete.push_back(block);
        } else {
          mfem::HypreParMatrix *block = new mfem::HypreParMatrix(*N_sparse_block);
          tmp_blocks(i, j) = block;
        }

      } else {
        MFEM_ABORT("Failed to cast operator blocks to mfem::HypreParMatrix");
      }
    }
  }
  mfem::HypreParMatrix *JJ(mfem::HypreParMatrixFromBlocks(tmp_blocks));
  Jacobian = JJ;
  for (auto ptr : blocks_to_delete) {
    delete ptr;
  }
  blocks_to_delete.clear();
  // std::unique_ptr<mfem::SparseMatrix> localJ(Add(1.0, M_->SpMat(), 0.0, M_->SpMat()));
  // localJ->Add(dt_, N_->GetLocalGradient(z));
  // Jacobian = M_->ParallelAssemble(localJ.get());
  // TODO(cci) adapt
  // std::unique_ptr<mfem::HypreParMatrix> Je(Jacobian->EliminateRowsCols(ess_tdof_list));

  return *Jacobian;
}

/**
 * @brief Destroy the Phase Field Reduced Operator:: Phase Field Reduced Operator object
 *
 */
PhaseFieldReducedOperator::~PhaseFieldReducedOperator() { delete Jacobian; }

/**
 * @file SteadyReducedOperator.hpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief Steady version of the linear system resulting from the NonLinear algorithm
 * @version 0.1
 * @date 2025-06-17
 *
 * @copyright Copyright (c) 2025
 *
 */
#include <memory>
#include <vector>

#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
#pragma once
/*
 *  Class SteadyPhaseFieldReducedOperator
 */
class SteadyPhaseFieldReducedOperator : public mfem::Operator {
 private:
  // RHS
  mfem::ParBlockNonlinearForm *RHS_;
  // Jacobian matrix
  mutable mfem::HypreParMatrix *Jacobian;

  // Time step
  double dt_;
  // Unknown
  const mfem::Vector *unk_;

  const std::vector<mfem::Array<int>> &ess_tdof_list;

 public:
  SteadyPhaseFieldReducedOperator(mfem::ParBlockNonlinearForm *RHS,
                                  const std::vector<mfem::Array<int>> &ess_tdof);

  /// Compute y = N(unk + dt*k) + M k
  void Mult(const mfem::Vector &k, mfem::Vector &y) const;

  /// Compute y = dt*grad_N(unk + dt*k) + M
  mfem::Operator &GetGradient(const mfem::Vector &k) const;
  ~SteadyPhaseFieldReducedOperator();
};

/**
 * @brief Construct a new Phase Field Reduced Operator:: Phase Field Reduced Operator object
 *
 * @param M
 * @param N
 */
SteadyPhaseFieldReducedOperator::SteadyPhaseFieldReducedOperator(
    mfem::ParBlockNonlinearForm *RHS, const std::vector<mfem::Array<int>> &ess_tdof)
    : Operator(RHS->Height()),
      RHS_(RHS),
      Jacobian(NULL),
      dt_(0.0),
      unk_(NULL),
      ess_tdof_list(ess_tdof) {}

/**
 * @brief  Compute y = N(unk + dt*k) + M k
 *
 * @param k
 * @param y
 */
void SteadyPhaseFieldReducedOperator::Mult(const mfem::Vector &k, mfem::Vector &y) const {
  this->RHS_->Mult(k, y);

  // TODO(cci) simplify BCs
  const mfem::Array<int> offsets = this->RHS_->GetBlockOffsets();
  const int fes_size = offsets.Size() - 1;
  auto sc_1 = 0;
  auto sc_2 = this->RHS_->Height() / fes_size;
  for (int i = 0; i < fes_size; ++i) {
    mfem::Vector y_i(y.GetData() + sc_1, sc_2);
    y_i.SetSubVector(ess_tdof_list[i], 0.0);
    sc_1 += sc_2;
  }
}

/**
 * @brief  Compute Jacobian
 *
 * @param k
 * @return mfem::Operator&
 */
mfem::Operator &SteadyPhaseFieldReducedOperator::GetGradient(const mfem::Vector &k) const {
  if (Jacobian != nullptr) {
    delete Jacobian;
  }
  const mfem::Array<int> offsets = this->RHS_->GetBlockOffsets();
  const int fes_size = offsets.Size() - 1;
  // Gets gradients of RHS_
  mfem::Operator &RHS_grad = this->RHS_->GetGradient(k);
  // Converts operators into BlockOperator
  mfem::BlockOperator *RHS_block_grad = dynamic_cast<mfem::BlockOperator *>(&RHS_grad);
  mfem::Array2D<mfem::HypreParMatrix *> tmp_blocks(fes_size, fes_size);
  std::vector<mfem::HypreParMatrix *> blocks_to_delete;

  for (int i = 0; i < fes_size; ++i) {
    for (int j = 0; j < fes_size; ++j) {
      mfem::Operator *RHS_block = &(RHS_block_grad->GetBlock(i, j));

      mfem::HypreParMatrix *RHS_sparse_block = dynamic_cast<mfem::HypreParMatrix *>(RHS_block);

      if (RHS_sparse_block) {
        mfem::HypreParMatrix *block = new mfem::HypreParMatrix(*RHS_sparse_block);
        // TODO(CCI) check_CI
        block->EliminateRowsCols(ess_tdof_list[i]);
        // CCI
        tmp_blocks(i, j) = block;

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
  return *Jacobian;
}

/**
 * @brief Destroy the Phase Field Reduced Operator:: Phase Field Reduced Operator object
 *
 */
SteadyPhaseFieldReducedOperator::~SteadyPhaseFieldReducedOperator() { delete Jacobian; }

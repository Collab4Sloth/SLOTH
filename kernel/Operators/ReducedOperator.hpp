/**
 * @file ReducedOperator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief UnSteady version of the linear system resulting from the NonLinear algorithm
 * @version 0.1
 * @date 2025-09-05
 *
 * Copyright CEA (C) 2025
 *
 * This file is part of SLOTH.
 *
 * SLOTH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SLOTH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include <vector>

#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
#pragma once
/*
 *  Class PhaseFieldReducedOperator
 */
class PhaseFieldReducedOperator : public mfem::Operator {
 private:
  // LHS
  mfem::ParBlockNonlinearForm* LHS_;
  // RHS
  mfem::ParBlockNonlinearForm* RHS_;
  // Jacobian matrix
  mutable std::unique_ptr<mfem::HypreParMatrix> Jacobian;
  // Time step
  double dt_;
  // Unknown
  const mfem::Vector* unk_;
  mutable mfem::Vector z;

  const std::vector<mfem::Array<int>>& ess_tdof_list;

  int fes_size_;

  mutable mfem::Array2D<const mfem::HypreParMatrix*> tmp_blocks_;
  mutable std::vector<std::unique_ptr<mfem::HypreParMatrix>> blocks_to_delete_;

 public:
  PhaseFieldReducedOperator(mfem::ParBlockNonlinearForm* M, mfem::ParBlockNonlinearForm* N,
                            const std::vector<mfem::Array<int>>& ess_tdof);

  /// Set current dt, unk values - needed to compute action and Jacobian.
  void SetParameters(double dt, const mfem::Vector* unk);

  /// Compute y = N(unk + dt*k) + M k
  void Mult(const mfem::Vector& k, mfem::Vector& y) const;

  /// Compute y = dt*grad_N(unk + dt*k) + M
  mfem::Operator& GetGradient(const mfem::Vector& k) const;
  ~PhaseFieldReducedOperator() = default;
};

/**
 * @brief Construct a new Phase Field Reduced Operator:: Phase Field Reduced Operator object
 *
 * @param LHS
 * @param RHS
 * @param ess_tdof
 */
PhaseFieldReducedOperator::PhaseFieldReducedOperator(mfem::ParBlockNonlinearForm* LHS,
                                                     mfem::ParBlockNonlinearForm* RHS,
                                                     const std::vector<mfem::Array<int>>& ess_tdof)
    : Operator(RHS->Height()),
      // : Operator(N->ParFESpace()->TrueVSize()),
      LHS_(LHS),
      RHS_(RHS),
      dt_(0.0),
      unk_(NULL),
      z(height),
      ess_tdof_list(ess_tdof) {
  const mfem::Array<int> offsets = this->RHS_->GetBlockOffsets();
  this->fes_size_ = offsets.Size() - 1;
  this->tmp_blocks_.SetSize(this->fes_size_, this->fes_size_);
  this->blocks_to_delete_.resize(this->fes_size_ * this->fes_size_);
}

/**
 * @brief  Set current dt, unk values - needed to compute action and Jacobian.
 *
 * @param dt
 * @param unk
 */
void PhaseFieldReducedOperator::SetParameters(double dt, const mfem::Vector* unk) {
  dt_ = dt;
  unk_ = unk;
}

/**
 * @brief  Compute y = N(unk + dt*k) + M k
 *
 * @param k
 * @param y
 */
void PhaseFieldReducedOperator::Mult(const mfem::Vector& k, mfem::Vector& y) const {
  add(*unk_, dt_, k, z);
  this->RHS_->Mult(z, y);
  this->LHS_->AddMult(k, y);

  // TODO(cci) simplify BCs
  auto sc_1 = 0;
  auto sc_2 = this->RHS_->Height() / this->fes_size_;
  for (int i = 0; i < this->fes_size_; ++i) {
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
mfem::Operator& PhaseFieldReducedOperator::GetGradient(const mfem::Vector& k) const {
  Jacobian.reset();
  UtilsForDebug::memory_checkpoint("avant getgradient Add");

  add(*unk_, dt_, k, z);
  // Gets gradients of RHS_ and LHS_
  UtilsForDebug::memory_checkpoint("avant getgradient LHS");
  mfem::BlockOperator& LHS_grad = this->LHS_->GetGradient(z);
  UtilsForDebug::memory_checkpoint("avant getgradient RHS");
  mfem::BlockOperator& RHS_grad = this->RHS_->GetGradient(z);

  UtilsForDebug::memory_checkpoint("avant assemblage RHS");
  for (int i = 0; i < this->fes_size_; ++i) {
    for (int j = 0; j < this->fes_size_; ++j) {
      const mfem::Operator& LHS_block = LHS_grad.GetBlock(i, j);
      const mfem::Operator& RHS_block = RHS_grad.GetBlock(i, j);

      const mfem::HypreParMatrix* LHS_sparse_block =
          dynamic_cast<const mfem::HypreParMatrix*>(&LHS_block);
      const mfem::HypreParMatrix* RHS_sparse_block =
          dynamic_cast<const mfem::HypreParMatrix*>(&RHS_block);

      if (!LHS_sparse_block || !RHS_sparse_block)
        MFEM_ABORT("Failed to cast operator blocks to mfem::HypreParMatrix");

      if (blocks_to_delete_[i * fes_size_ + j]) {
        blocks_to_delete_[i * fes_size_ + j].reset();
      }
      blocks_to_delete_[i * fes_size_ + j] =
          std::make_unique<mfem::HypreParMatrix>(*LHS_sparse_block);

      blocks_to_delete_[i * fes_size_ + j]->Add(dt_, *RHS_sparse_block);

      std::unique_ptr<mfem::HypreParMatrix> bb(
          blocks_to_delete_[i * fes_size_ + j]->EliminateRowsCols(ess_tdof_list[i]));

      tmp_blocks_(i, j) = blocks_to_delete_[i * fes_size_ + j].get();
    }
  }
  Jacobian.reset(mfem::HypreParMatrixFromBlocks(tmp_blocks_));
  UtilsForDebug::memory_checkpoint("avant return jacobian ");

  return *Jacobian;
}

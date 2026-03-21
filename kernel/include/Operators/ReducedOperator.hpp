/**
 * @file ReducedOperator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief UnSteady version of the linear system resulting from the NonLinear algorithm
 * @version 0.1
 * @date 2025-09-05
 *
 * @copyright CEA (C) 2025
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
#pragma once
#include <memory>
#include <vector>

#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

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

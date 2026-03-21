
/**
 * @file ISolverBase.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief  Base class for iterative solvers
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
#include <string>

#include "Options/Options.hpp"
#include "Parameters/Parameters.hpp"
#include "Solvers/SolverBase.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

class SolverBICGSTAB : public SolverBase<mfem::BiCGSTABSolver, IterativeSolverType> {
 public:
  SolverBICGSTAB() {}
  std::shared_ptr<mfem::BiCGSTABSolver> create_solver(const Parameters& params) override;

  ~SolverBICGSTAB() {}
};

////////////////////////////////////
////////////////////////////////////

class SolverMINRES : public SolverBase<mfem::MINRESSolver, IterativeSolverType> {
 public:
  SolverMINRES() {}
  std::shared_ptr<mfem::MINRESSolver> create_solver(const Parameters& params) override;

  ~SolverMINRES() {}
};

////////////////////////////////////
////////////////////////////////////

class SolverCG : public SolverBase<mfem::CGSolver, IterativeSolverType> {
 public:
  SolverCG() {}
  std::shared_ptr<mfem::CGSolver> create_solver(const Parameters& params) override;

  ~SolverCG() {}
};
////////////////////////////////////
////////////////////////////////////

class SolverGMRES : public SolverBase<mfem::GMRESSolver, IterativeSolverType> {
 public:
  SolverGMRES() {}
  std::shared_ptr<mfem::GMRESSolver> create_solver(const Parameters& params) override;

  ~SolverGMRES() {}
};

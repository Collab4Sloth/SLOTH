
/**
 * @file NLSolverBase.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief ase class for non linear solvers
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
#include <memory>
#include <string>

#include "Options/Options.hpp"
#include "Solvers/SolverBase.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

class NLSolverBase : SolverBase<mfem::NewtonSolver, NLSolverType> {
 public:
  NLSolverBase();
  std::shared_ptr<mfem::NewtonSolver> create_solver(NLSolverType SOLVER,
                                                    const Parameters& params) override;

  ~NLSolverBase();
};

/**
 * @brief Construct a new NLSolverBase::NLSolverBase object
 *
 */
NLSolverBase::NLSolverBase() {}

/**
 * @brief Create non linear solver based of the NLSolverType and a list of Parameters
 *
 * @param NLSOLVER
 * @param params
 * @return std::shared_ptr<mfem::NewtonSolver>
 */
std::shared_ptr<mfem::NewtonSolver> NLSolverBase::create_solver(NLSolverType NLSOLVER,
                                                                const Parameters& params) {
  // TODO(cci): mettre un check sur le param de base name
  // TODO(cc) : mettre un getinfo pour la doc
  this->solver_description_ = params.get_param_value<std::string>("description");
  switch (NLSOLVER) {
    case NLSolverType::NEWTON: {
      const int print_level =
          params.get_param_value_or_default<int>("print_level", NewtonDefaultConstant::print_level);
      const bool iterative_mode = params.get_param_value_or_default<bool>(
          "iterative_mode", NewtonDefaultConstant::iterative_mode);
      const int n_iter_max =
          params.get_param_value_or_default<int>("iter_max", NewtonDefaultConstant::iter_max);
      const double n_rel_tol =
          params.get_param_value_or_default<double>("rel_tol", NewtonDefaultConstant::rel_tol);
      const double n_abs_tol =
          params.get_param_value_or_default<double>("abs_tol", NewtonDefaultConstant::abs_tol);

      auto ss = std::make_shared<mfem::NewtonSolver>(MPI_COMM_WORLD);
      // mfem::KINSolver* kinsolver = new mfem::KINSolver(MPI_COMM_WORLD, 1);
      // kinsolver->SetJFNK(true);
      // kinsolver->SetLSMaxter(100);
      // ss = kinsolver;

      ss->SetPrintLevel(print_level);
      ss->iterative_mode = iterative_mode;
      ss->SetMaxIter(n_iter_max);
      ss->SetRelTol(n_rel_tol);
      ss->SetAbsTol(n_abs_tol);
      // kinsolver->SetMaxSetupCalls(4);

      return ss;
    }
    default:
      mfem::mfem_error("Solver::create_precond: SMOOTHER, UMFPACK are available");
      break;
  }
}
/**
 * @brief Destroy the NLSolverBase::NLSolverBase object
 *
 */
NLSolverBase::~NLSolverBase() {}

/**
 * @file LSolver.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Class used to manage Linear Solver  objet
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

#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Solvers/SlothSolver.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

/**
 * @brief Generic class for defining linear solvers with or without preconditionners
 *
 */
class LSolver {
 private:
  VSharedMFEMSolver variant_solver_;
  std::shared_ptr<SlothSolver> ss;

  VSharedMFEMSolver variant_precond_;
  std::shared_ptr<SlothSolver> pp;

 public:
  LSolver(VSolverType SOLVER, const Parameters& s_params, mfem::Operator& ope);
  LSolver(VSolverType SOLVER, const Parameters& s_params, VSolverType PRECOND,
          const Parameters& p_params, mfem::Operator& ope);
  VSharedMFEMSolver get_solver();

  ~LSolver();
};

/**
 * @brief Construct a new LSolver::LSolver object with a preconditionner
 *
 * @param SOLVER
 * @param s_params
 * @param PRECOND
 * @param p_params
 * @param ope
 */
LSolver::LSolver(VSolverType SOLVER, const Parameters& s_params, VSolverType PRECOND,
                 const Parameters& p_params, mfem::Operator& ope) {
  SlothInfo::debug("LSolver::LSolver start");
  ss = std::make_shared<SlothSolver>(SOLVER, s_params);
  this->variant_solver_ = ss->get_value();

  pp = std::make_shared<SlothSolver>(PRECOND, p_params);
  this->variant_precond_ = pp->get_value();

  SetPrecondSolver func_prec = {ope};
  std::visit(func_prec, this->variant_solver_, this->variant_precond_);
  SlothInfo::debug("LSolver::LSolver end");
}

/**
 * @brief Construct a new LSolver::LSolver object without precondnitioner
 *
 * @param SOLVER
 * @param s_params
 * @param ope
 */
LSolver::LSolver(VSolverType SOLVER, const Parameters& s_params, mfem::Operator& ope) {
  SlothInfo::debug("LSolver::LSolver start");

  ss = std::make_shared<SlothSolver>(SOLVER, s_params);
  this->variant_solver_ = ss->get_value();

  this->variant_precond_ = std::make_shared<std::monostate>();

  SetPrecondSolver func_prec = {ope};
  std::visit(func_prec, this->variant_solver_, this->variant_precond_);

  SlothInfo::debug("LSolver::LSolver end");
}

/**
 * @brief Return the solver
 *
 * @return VSharedMFEMSolver
 */
VSharedMFEMSolver LSolver::get_solver() { return this->variant_solver_; }

/**
 * @brief Destroy the LSolver::LSolver object
 *
 */
LSolver::~LSolver() {}

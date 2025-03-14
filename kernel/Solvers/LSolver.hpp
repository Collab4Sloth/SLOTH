/**
 * @file LSolver.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Class used to manage Linear Solver  objet
 * @version 0.1
 * @date 2024-05-25
 *
 * Copyright CEA (c) 2024
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
  VSharedMFEMSolver get_solver() { return this->variant_solver_; }
  ~LSolver() = default;
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
DEBILE_INLINE LSolver::LSolver(VSolverType SOLVER, const Parameters& s_params, VSolverType PRECOND,
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
DEBILE_INLINE LSolver::LSolver(VSolverType SOLVER, const Parameters& s_params, mfem::Operator& ope) {
  SlothInfo::debug("LSolver::LSolver start");

  ss = std::make_shared<SlothSolver>(SOLVER, s_params);
  this->variant_solver_ = ss->get_value();

  this->variant_precond_ = std::make_shared<std::monostate>();

  SetPrecondSolver func_prec = {ope};
  std::visit(func_prec, this->variant_solver_, this->variant_precond_);

  SlothInfo::debug("LSolver::LSolver end");
}

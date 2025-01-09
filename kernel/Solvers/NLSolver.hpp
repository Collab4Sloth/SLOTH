/**
 * @file NLSolver.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Class used to manage a NonLinear Solver( Newton type)
 * @version 0.1
 * @date 2024-07-24
 *
 * Copyright CEA (c) 2024
 *
 */
#include <memory>
#include <type_traits>
#include <variant>

#include "Options/Options.hpp"
#include "Solvers/DSolverBase.hpp"
#include "Solvers/IPrecondBase.hpp"
#include "Solvers/ISolverBase.hpp"
#include "Solvers/NLSolverBase.hpp"
#include "Solvers/SlothSolver.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

class NLSolver {
 private:
  VSharedMFEMSolver variant_solver_;
  std::shared_ptr<SlothSolver> ss;

  VSharedMFEMSolver variant_precond_;
  std::shared_ptr<SlothSolver> pp;

  NLSolverBase NLSolverBase_;
  // TODO(cci) on part sur newton mais on pourrait compléter la liste
  std::shared_ptr<mfem::NewtonSolver> nl_solver_;

 public:
  NLSolver(NLSolverType NLSOLVER, const Parameters& nl_params, VSolverType SOLVER,
           const Parameters& s_params, VSolverType PRECOND, const Parameters& p_params,
           mfem::Operator& ope);
  NLSolver(NLSolverType NLSOLVER, const Parameters& nl_params, VSolverType SOLVER,
           const Parameters& s_params, mfem::Operator& ope);
  std::shared_ptr<mfem::NewtonSolver> get_nl_solver();

  ~NLSolver();
};

/**
 * @brief Construct a new NLSolver::NLSolver object
 *
 * @param NLSOLVER
 * @param nl_params
 * @param SOLVER
 * @param s_params
 * @param PRECOND
 * @param p_params
 * @param ope
 */
NLSolver::NLSolver(NLSolverType NLSOLVER, const Parameters& nl_params, VSolverType SOLVER,
                   const Parameters& s_params, VSolverType PRECOND, const Parameters& p_params,
                   mfem::Operator& ope)
    : nl_solver_(NLSolverBase_.create_solver(NLSOLVER, nl_params)) {
  SlothInfo::debug("NLSolver::NLSolver start");
  ss = std::make_shared<SlothSolver>(SOLVER, s_params);
  this->variant_solver_ = ss->get_value();

  pp = std::make_shared<SlothSolver>(PRECOND, p_params);
  this->variant_precond_ = pp->get_value();

  SetPrecondNLSolver func_prec = {ope, nl_solver_};
  std::visit(func_prec, this->variant_solver_, this->variant_precond_);
}

/**
 * @brief Construct a new NLSolver::NLSolver object
 *
 * @param NLSOLVER
 * @param nl_params
 * @param SOLVER
 * @param s_params
 * @param ope
 */
NLSolver::NLSolver(NLSolverType NLSOLVER, const Parameters& nl_params, VSolverType SOLVER,
                   const Parameters& s_params, mfem::Operator& ope)
    : nl_solver_(NLSolverBase_.create_solver(NLSOLVER, nl_params)) {
  SlothInfo::debug("NLSolver::NLSolver start");
  ss = std::make_shared<SlothSolver>(SOLVER, s_params);
  this->variant_solver_ = ss->get_value();

  this->variant_precond_ = std::make_shared<std::monostate>();

  SetPrecondNLSolver func_prec = {ope, nl_solver_};
  std::visit(func_prec, this->variant_solver_, this->variant_precond_);
}

/**
 * @brief Return the Non Linear Solver
 *
 * @return std::shared_ptr<mfem::NewtonSolver>
 */
std::shared_ptr<mfem::NewtonSolver> NLSolver::get_nl_solver() { return this->nl_solver_; }

/**
 * @brief Destroy the Solver::Solver object
 *
 */
NLSolver::~NLSolver() {}

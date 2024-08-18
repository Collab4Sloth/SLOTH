
/**
 * @file NLSolverBase.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief  Base class for non linear solvers
 * @version 0.1
 * @date 2024-07-25
 *
 * Copyright CEA (c) 2024
 *
 */
#include <memory>
#include <string>

#include "Solvers/SolverBase.hpp"
#include "Utils/PhaseFieldConstants.hpp"
#include "Utils/PhaseFieldOptions.hpp"
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
      ss->SetPrintLevel(print_level);
      ss->iterative_mode = iterative_mode;
      ss->SetMaxIter(n_iter_max);
      ss->SetRelTol(n_rel_tol);
      ss->SetAbsTol(n_abs_tol);
      return ss;
    }
    default:
      throw std::runtime_error("Solver::create_precond: SMOOTHER, UMFPACK are available");
      break;
  }
}
/**
 * @brief Destroy the NLSolverBase::NLSolverBase object
 *
 */
NLSolverBase::~NLSolverBase() {}

/*
 * UtilsForSolvers.hpp
 *
 *  Created on: 20/3/2023
 *      Author: ci230846
 */
#include "mfem.hpp"

#ifndef UTILSFORSOLVERS_H_
#define UTILSFORSOLVERS_H_
// TODO(ci230846) : ajouter une gestion des solvers avec des enums.

/**
 * @brief Useful methods for managing solvers
 *
 */
class UtilsForSolvers {
 public:
  void BuildSolver(mfem::IterativeSolver& solver, mfem::Solver& solv_method, mfem::Operator& ope);
  void SetSolverParameters(mfem::IterativeSolver& solver, const int& _print_level,
                           const bool _iterative_mode, const int& _n_iter_max,
                           const double& _n_rel_tol, const double& _n_abs_tol);
};

/**
 * @brief Build solver depending on given preconditionner and operator
 *
 * @param solver IterativeSolver
 * @param preconditionner Preconditionner
 * @param ope Operator
 */
void UtilsForSolvers::BuildSolver(mfem::IterativeSolver& solver, mfem::Solver& preconditionner,
                                  mfem::Operator& ope) {
  solver.SetPreconditioner(preconditionner);
  solver.SetOperator(ope);
}
/**
 * @brief Set all parameters used by the solver
 *
 * @param solver IterativeSolver
 * @param _print_level printing level
 * @param _iterative_mode flag to activate the iterative mode (initialization of the Iterative
 * Solver, False by default)
 * @param _n_iter_max maximum number of iterations used by the IterativeSolver
 * @param _n_rel_tol relative tolerance used by the IterativeSolver convergence test
 * @param _n_abs_tol absolute tolerance used by the IterativeSolver convergence test
 */
void UtilsForSolvers::SetSolverParameters(mfem::IterativeSolver& solver, const int& _print_level,
                                          const bool _iterative_mode, const int& _n_iter_max,
                                          const double& _n_rel_tol, const double& _n_abs_tol) {
  solver.SetPrintLevel(_print_level);
  solver.iterative_mode = _iterative_mode;
  solver.SetMaxIter(_n_iter_max);
  solver.SetRelTol(_n_rel_tol);
  solver.SetAbsTol(_n_abs_tol);
}
#endif /* UTILSFORSOLVERS_H_ */


/**
 * @file ISolverBase.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief  Base class for iterative solvers
 * @version 0.1
 * @date 2024-07-25
 *
 * Copyright CEA (c) 2024
 *
 */
#include <memory>
#include <string>

#include "Options/Options.hpp"
#include "Solvers/SolverBase.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

class SolverBICGSTAB : public SolverBase<mfem::BiCGSTABSolver, IterativeSolverType> {
 public:
  SolverBICGSTAB() {}
  std::shared_ptr<mfem::BiCGSTABSolver> create_solver(IterativeSolverType SOLVER,
                                                      const Parameters& params) override;

  ~SolverBICGSTAB() {}
};

/**
 * @brief Create a iterative solver based of the SolverType and a list of Parameters
 *
 * @param SOLVER
 * @param params
 * @return std::shared_ptr<mfem::IterativeSolver>
 */
std::shared_ptr<mfem::BiCGSTABSolver> SolverBICGSTAB::create_solver(IterativeSolverType SOLVER,
                                                                    const Parameters& params) {
  // TODO(cci): mettre un check sur le param de base name
  // TODO(cc) : mettre un getinfo pour la doc
  this->solver_description_ = params.get_param_value<std::string>("description");
  SlothInfo::debug(" Create ", this->get_description());

  const int print_level =
      params.get_param_value_or_default<int>("print_level", BICGSTABDefaultConstant::print_level);
  const bool iterative_mode = params.get_param_value_or_default<bool>(
      "iterative_mode", BICGSTABDefaultConstant::iterative_mode);
  const int n_iter_max =
      params.get_param_value_or_default<int>("iter_max", BICGSTABDefaultConstant::iter_max);
  const double n_rel_tol =
      params.get_param_value_or_default<double>("rel_tol", BICGSTABDefaultConstant::rel_tol);
  const double n_abs_tol =
      params.get_param_value_or_default<double>("abs_tol", BICGSTABDefaultConstant::abs_tol);

  auto ss = std::make_shared<mfem::BiCGSTABSolver>(MPI_COMM_WORLD);
  ss->SetPrintLevel(print_level);
  ss->iterative_mode = iterative_mode;
  ss->SetMaxIter(n_iter_max);
  ss->SetRelTol(n_rel_tol);
  ss->SetAbsTol(n_abs_tol);
  return ss;
}

////////////////////////////////////
////////////////////////////////////

class SolverMINRES : public SolverBase<mfem::MINRESSolver, IterativeSolverType> {
 public:
  SolverMINRES() {}
  std::shared_ptr<mfem::MINRESSolver> create_solver(IterativeSolverType SOLVER,
                                                    const Parameters& params) override;

  ~SolverMINRES() {}
};

std::shared_ptr<mfem::MINRESSolver> SolverMINRES::create_solver(IterativeSolverType SOLVER,
                                                                const Parameters& params) {
  // TODO(cci): mettre un check sur le param de base name
  // TODO(cc) : mettre un getinfo pour la doc
  this->solver_description_ = params.get_param_value<std::string>("description");
  SlothInfo::debug(" Create ", this->get_description());

  const int print_level =
      params.get_param_value_or_default<int>("print_level", MINRESDefaultConstant::print_level);
  const bool iterative_mode = params.get_param_value_or_default<bool>(
      "iterative_mode", MINRESDefaultConstant::iterative_mode);
  const int n_iter_max =
      params.get_param_value_or_default<int>("iter_max", MINRESDefaultConstant::iter_max);
  const double n_rel_tol =
      params.get_param_value_or_default<double>("rel_tol", MINRESDefaultConstant::rel_tol);
  const double n_abs_tol =
      params.get_param_value_or_default<double>("abs_tol", MINRESDefaultConstant::abs_tol);

  auto ss = std::make_shared<mfem::MINRESSolver>(MPI_COMM_WORLD);
  ss->SetPrintLevel(print_level);
  ss->iterative_mode = iterative_mode;
  ss->SetMaxIter(n_iter_max);
  ss->SetRelTol(n_rel_tol);
  ss->SetAbsTol(n_abs_tol);
  return ss;
}

////////////////////////////////////
////////////////////////////////////

class SolverCG : public SolverBase<mfem::CGSolver, IterativeSolverType> {
 public:
  SolverCG() {}
  std::shared_ptr<mfem::CGSolver> create_solver(IterativeSolverType SOLVER,
                                                const Parameters& params) override;

  ~SolverCG() {}
};
std::shared_ptr<mfem::CGSolver> SolverCG::create_solver(IterativeSolverType SOLVER,
                                                        const Parameters& params) {
  // TODO(cci): mettre un check sur le param de base name
  // TODO(cc) : mettre un getinfo pour la doc
  this->solver_description_ = params.get_param_value<std::string>("description");
  SlothInfo::debug(" Create ", this->get_description());

  const int print_level =
      params.get_param_value_or_default<int>("print_level", CGDefaultConstant::print_level);
  const bool iterative_mode =
      params.get_param_value_or_default<bool>("iterative_mode", CGDefaultConstant::iterative_mode);
  const int n_iter_max =
      params.get_param_value_or_default<int>("iter_max", CGDefaultConstant::iter_max);
  const double n_rel_tol =
      params.get_param_value_or_default<double>("rel_tol", CGDefaultConstant::rel_tol);
  const double n_abs_tol =
      params.get_param_value_or_default<double>("abs_tol", CGDefaultConstant::abs_tol);

  auto ss = std::make_shared<mfem::CGSolver>(MPI_COMM_WORLD);
  ss->SetPrintLevel(print_level);
  ss->iterative_mode = iterative_mode;
  ss->SetMaxIter(n_iter_max);
  ss->SetRelTol(n_rel_tol);
  ss->SetAbsTol(n_abs_tol);
  return ss;
}

////////////////////////////////////
////////////////////////////////////

class SolverGMRES : public SolverBase<mfem::GMRESSolver, IterativeSolverType> {
 public:
  SolverGMRES() {}
  std::shared_ptr<mfem::GMRESSolver> create_solver(IterativeSolverType SOLVER,
                                                   const Parameters& params) override;

  ~SolverGMRES() {}
};

std::shared_ptr<mfem::GMRESSolver> SolverGMRES::create_solver(IterativeSolverType SOLVER,
                                                              const Parameters& params) {
  // TODO(cci): mettre un check sur le param de base name
  // TODO(cc) : mettre un getinfo pour la doc
  this->solver_description_ = params.get_param_value<std::string>("description");
  SlothInfo::debug(" Create ", this->get_description());

  const int print_level =
      params.get_param_value_or_default<int>("print_level", GMRESDefaultConstant::print_level);
  const bool iterative_mode = params.get_param_value_or_default<bool>(
      "iterative_mode", GMRESDefaultConstant::iterative_mode);
  const int n_iter_max =
      params.get_param_value_or_default<int>("iter_max", GMRESDefaultConstant::iter_max);
  const double n_rel_tol =
      params.get_param_value_or_default<double>("rel_tol", GMRESDefaultConstant::rel_tol);
  const double n_abs_tol =
      params.get_param_value_or_default<double>("abs_tol", GMRESDefaultConstant::abs_tol);

  auto ss = std::make_shared<mfem::GMRESSolver>(MPI_COMM_WORLD);
  ss->SetPrintLevel(print_level);
  ss->iterative_mode = iterative_mode;
  ss->SetMaxIter(n_iter_max);
  ss->SetRelTol(n_rel_tol);
  ss->SetAbsTol(n_abs_tol);
  return ss;
}

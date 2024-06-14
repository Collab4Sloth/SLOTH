/**
 * @file Solver.cpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Class to define a Solver objet
 * @version 0.1
 * @date 2024-05-25
 *
 * Copyright CEA (c) 2024
 *
 */

#include "Utils/PhaseFieldConstants.hpp"
#include "Utils/PhaseFieldOptions.hpp"
#include "mfem.hpp"

#pragma once

class Solver {
 private:
  std::shared_ptr<mfem::IterativeSolver> create_solver(SolverType SOLVER);
  std::shared_ptr<mfem::Solver> create_precond(PreconditionerType PRECON);

  std::shared_ptr<mfem::IterativeSolver> solver_;
  std::shared_ptr<mfem::Solver> precond_;
  int print_level_;
  bool iterative_mode_;
  int n_iter_max_;
  double n_rel_tol_;
  double n_abs_tol_;
  int prec_print_level_;
  bool prec_iterative_mode_;
  int prec_n_iter_max_;
  double prec_n_rel_tol_;
  double prec_n_abs_tol_;
  bool get_default_{false};
  bool prec_get_default_{false};

  void set_solver_parameters();

 public:
  Solver(SolverType SOLVER, PreconditionerType PRECOND, mfem::Operator& ope);

  Solver(SolverType SOLVER, const int& print_level, const bool iterative_mode,
         const int& n_iter_max, const double& n_rel_tol, const double& n_abs_tol,
         PreconditionerType PRECOND, mfem::Operator& ope);

  Solver(SolverType SOLVER, PreconditionerType PRECOND, const int& prec_print_level,
         const bool prec_iterative_mode, const int& prec_n_iter_max, const double& prec_n_rel_tol,
         const double& prec_n_abs_tol, mfem::Operator& ope);

  Solver(SolverType SOLVER, const int& print_level, const bool iterative_mode,
         const int& n_iter_max, const double& n_rel_tol, const double& n_abs_tol,
         PreconditionerType PRECOND, const int& prec_print_level, const bool prec_iterative_mode,
         const int& prec_n_iter_max, const double& prec_n_rel_tol, const double& prec_n_abs_tol,
         mfem::Operator& ope);
  std::shared_ptr<mfem::Solver> get_precond();
  std::shared_ptr<mfem::IterativeSolver> get_solver();

  ~Solver();
};

/**
 * @brief  Build the solver depending on SolverType value
 *
 * @param SOLVER
 * @return std::shared_ptr<mfem::IterativeSolver>
 */
std::shared_ptr<mfem::IterativeSolver> Solver::create_solver(SolverType SOLVER) {
  switch (SOLVER) {
    case SolverType::NEWTON: {
      if (this->get_default_) {
        this->print_level_ = NewtonDefaultConstant::print_level;
        this->iterative_mode_ = NewtonDefaultConstant::iterative_mode;
        this->n_iter_max_ = NewtonDefaultConstant::iter_max;
        this->n_rel_tol_ = NewtonDefaultConstant::rel_tol;
        this->n_abs_tol_ = NewtonDefaultConstant::abs_tol;
      }
      return std::make_shared<mfem::NewtonSolver>();
    }
    case SolverType::BICGSTAB: {
      if (this->get_default_) {
        this->print_level_ = BICGSTABDefaultConstant::print_level;
        this->iterative_mode_ = BICGSTABDefaultConstant::iterative_mode;
        this->n_iter_max_ = BICGSTABDefaultConstant::iter_max;
        this->n_rel_tol_ = BICGSTABDefaultConstant::rel_tol;
        this->n_abs_tol_ = BICGSTABDefaultConstant::abs_tol;
      }
      return std::make_shared<mfem::BiCGSTABSolver>();
    }
    case SolverType::CG: {
      if (this->get_default_) {
        this->print_level_ = CGDefaultConstant::print_level;
        this->iterative_mode_ = CGDefaultConstant::iterative_mode;
        this->n_iter_max_ = CGDefaultConstant::iter_max;
        this->n_rel_tol_ = CGDefaultConstant::rel_tol;
        this->n_abs_tol_ = CGDefaultConstant::abs_tol;
      }
      return std::make_shared<mfem::CGSolver>();
    }
    case SolverType::GMRES: {
      if (this->get_default_) {
        this->print_level_ = GMRESDefaultConstant::print_level;
        this->iterative_mode_ = GMRESDefaultConstant::iterative_mode;
        this->n_iter_max_ = GMRESDefaultConstant::iter_max;
        this->n_rel_tol_ = GMRESDefaultConstant::rel_tol;
        this->n_abs_tol_ = GMRESDefaultConstant::abs_tol;
      }
      return std::make_shared<mfem::GMRESSolver>();
    }
    default:
      throw std::runtime_error("Solver::create_solver: NEWTON are available");
      break;
  }
}

/**
 * @brief Build the preconditionner depending on PreconditionerType value
 *
 * @param PRECOND
 * @return std::shared_ptr<mfem::Solver>
 */
std::shared_ptr<mfem::Solver> Solver::create_precond(PreconditionerType PRECOND) {
  switch (PRECOND) {
    case PreconditionerType::SMOOTHER: {
      if (this->prec_get_default_) {
        this->prec_print_level_ = MassDefaultConstant::print_level;
        this->prec_iterative_mode_ = MassDefaultConstant::iterative_mode;
        this->prec_n_iter_max_ = MassDefaultConstant::iter_max;
        this->prec_n_rel_tol_ = MassDefaultConstant::rel_tol;
        this->prec_n_abs_tol_ = MassDefaultConstant::abs_tol;
      }
      return std::make_shared<mfem::DSmoother>();
    }
    case PreconditionerType::UMFPACK: {
      if (this->prec_get_default_) {
        this->prec_print_level_ = MassDefaultConstant::print_level;
        this->prec_iterative_mode_ = MassDefaultConstant::iterative_mode;
        this->prec_n_iter_max_ = MassDefaultConstant::iter_max;
        this->prec_n_rel_tol_ = MassDefaultConstant::rel_tol;
        this->prec_n_abs_tol_ = MassDefaultConstant::abs_tol;
      }
      return std::make_shared<mfem::UMFPackSolver>();
    }
    default:
      throw std::runtime_error("Solver::create_precond: SMOOTHER, UMFPACK are available");
      break;
  }
}

/**
 * @brief Construct a new Solver:: Solver object with default precondition's and solvers's
 * parameters
 *
 * @param SOLVER
 * @param PRECOND
 * @param ope
 */
Solver::Solver(SolverType SOLVER, PreconditionerType PRECOND, mfem::Operator& ope)
    : get_default_(true), prec_get_default_(true) {
  this->solver_ = this->create_solver(SOLVER);
  this->precond_ = this->create_precond(PRECOND);
  this->set_solver_parameters();

  auto prec = this->precond_.get();
  this->solver_->SetPreconditioner(*prec);
  this->solver_->SetOperator(ope);
}

/**
 * @brief Construct a new Solver:: Solver object with default preconditioner's parameters and given
 * solver's parameters
 * @param solver IterativeSolver
 * @param preconditionner Preconditionner
 * @param ope Operator
 */
Solver::Solver(SolverType SOLVER, const int& print_level, const bool iterative_mode,
               const int& n_iter_max, const double& n_rel_tol, const double& n_abs_tol,
               PreconditionerType PRECOND, mfem::Operator& ope)
    : prec_get_default_(true),
      print_level_(print_level),
      iterative_mode_(iterative_mode),
      n_iter_max_(n_iter_max),
      n_rel_tol_(n_rel_tol),
      n_abs_tol_(n_abs_tol) {
  this->solver_ = this->create_solver(SOLVER);
  this->precond_ = this->create_precond(PRECOND);
  this->set_solver_parameters();

  auto prec = this->precond_.get();
  this->solver_->SetPreconditioner(*prec);
  this->solver_->SetOperator(ope);
}

/**
 * @brief Construct a new Solver:: Solver object with default solver's parameters and given
 * preconditioner's parameters
 *
 * @param SOLVER
 * @param PRECOND
 * @param prec_print_level
 * @param prec_iterative_mode
 * @param prec_n_iter_max
 * @param prec_n_rel_tol
 * @param prec_n_abs_tol
 * @param ope
 */
Solver::Solver(SolverType SOLVER, PreconditionerType PRECOND, const int& prec_print_level,
               const bool prec_iterative_mode, const int& prec_n_iter_max,
               const double& prec_n_rel_tol, const double& prec_n_abs_tol, mfem::Operator& ope)
    : get_default_(true),
      prec_print_level_(prec_print_level),
      prec_iterative_mode_(prec_iterative_mode),
      prec_n_iter_max_(prec_n_iter_max),
      prec_n_rel_tol_(prec_n_rel_tol),
      prec_n_abs_tol_(prec_n_abs_tol) {
  this->solver_ = this->create_solver(SOLVER);
  this->precond_ = this->create_precond(PRECOND);
  this->set_solver_parameters();

  auto prec = this->precond_.get();
  this->solver_->SetPreconditioner(*prec);
  this->solver_->SetOperator(ope);
}

/**
 * @brief Construct a new Solver:: Solver object with given preconditioner's and solver's parameters
 *
 * @param SOLVER
 * @param print_level
 * @param iterative_mode
 * @param n_iter_max
 * @param n_rel_tol
 * @param n_abs_tol
 * @param PRECOND
 * @param prec_print_level
 * @param prec_iterative_mode
 * @param prec_n_iter_max
 * @param prec_n_rel_tol
 * @param prec_n_abs_tol
 * @param ope
 */
Solver::Solver(SolverType SOLVER, const int& print_level, const bool iterative_mode,
               const int& n_iter_max, const double& n_rel_tol, const double& n_abs_tol,
               PreconditionerType PRECOND, const int& prec_print_level,
               const bool prec_iterative_mode, const int& prec_n_iter_max,
               const double& prec_n_rel_tol, const double& prec_n_abs_tol, mfem::Operator& ope)
    : print_level_(print_level),
      iterative_mode_(iterative_mode),
      n_iter_max_(n_iter_max),
      n_rel_tol_(n_rel_tol),
      n_abs_tol_(n_abs_tol),
      prec_print_level_(prec_print_level),
      prec_iterative_mode_(prec_iterative_mode),
      prec_n_iter_max_(prec_n_iter_max),
      prec_n_rel_tol_(prec_n_rel_tol),
      prec_n_abs_tol_(prec_n_abs_tol) {
  this->solver_ = this->create_solver(SOLVER);
  this->precond_ = this->create_precond(PRECOND);
  this->set_solver_parameters();

  auto prec = this->precond_.get();
  this->solver_->SetPreconditioner(*prec);
  this->solver_->SetOperator(ope);
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
void Solver::set_solver_parameters() {
  this->solver_->SetPrintLevel(this->print_level_);
  this->solver_->iterative_mode = this->iterative_mode_;
  this->solver_->SetMaxIter(this->n_iter_max_);
  this->solver_->SetRelTol(this->n_rel_tol_);
  this->solver_->SetAbsTol(this->n_abs_tol_);
}

/**
 * @brief Return the preconditionner
 *
 * @return std::shared_ptr<mfem::Solver>
 */
std::shared_ptr<mfem::Solver> Solver::get_precond() { return this->precond_; }

/**
 * @brief Return the solver
 *
 * @return std::shared_ptr<mfem::IterativeSolver>
 */
std::shared_ptr<mfem::IterativeSolver> Solver::get_solver() { return this->solver_; }

/**
 * @brief Destroy the Solver:: Solver object
 *
 */
Solver::~Solver() {}
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

#include "Solvers/DSolverBase.hpp"
#include "Solvers/IPrecondBase.hpp"
#include "Solvers/ISolverBase.hpp"
#include "Solvers/NLSolverBase.hpp"
#include "Solvers/SlothSolver.hpp"
#include "Utils/PhaseFieldConstants.hpp"
#include "Utils/PhaseFieldOptions.hpp"
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

  std::visit(
      [this](auto&& arg) {
        using TT = std::decay_t<decltype(arg)>;
        if constexpr (is_in_variant_v<TT, VIterativeSolver>) {
          std::visit(
              [&](auto&& prec) {
                using PP = std::decay_t<decltype(prec)>;
                if constexpr (!std::is_same_v<PP, std::shared_ptr<std::monostate>>) {
                  MFEM_VERIFY((is_in_variant_v<PP, VIterativePrecond>),
                              "NLSolver:: IterativeSolver objects  can only be associated with an "
                              "IterativePreconditionner objects");
                  if constexpr (is_in_variant_v<PP, VIterativePrecond>) {
                    SlothInfo::debug("NLSolver::NLSolver setting preconditionner");

                    std::shared_ptr<mfem::Solver> aaPrec =
                        std::dynamic_pointer_cast<mfem::Solver>(prec);
                    arg->SetPreconditioner(*aaPrec);
                  }
                }
              },
              this->variant_precond_);
        }
        if constexpr (is_in_variant_v<TT, VHypreSolver>) {
          std::visit(
              [&](auto&& prec) {
                using PP = std::decay_t<decltype(prec)>;

                if constexpr (!std::is_same_v<PP, std::shared_ptr<std::monostate>>) {
                  MFEM_VERIFY((is_in_variant_v<PP, VHyprePrecond>),
                              "NLSolver:: HypreSolver objects  can only be associated with an "
                              "HyprePreconditionner objects");

                  if constexpr (is_in_variant_v<PP, VHyprePrecond>) {
                    SlothInfo::debug(
                        "NLSolver::NLSolver setting preconditionner (not hypre smoother)");

                    std::shared_ptr<mfem::HypreSolver> aaPrec =
                        std::dynamic_pointer_cast<mfem::HypreSolver>(prec);
                    arg->SetPreconditioner(*aaPrec);
                  }
                }
              },
              this->variant_precond_);
        }
      },
      this->variant_solver_);

  SlothInfo::debug("NLSolver::NLSolver setting operator ");
  this->nl_solver_->SetOperator(ope);

  std::visit(
      [this](auto&& arg) {
        using TT = std::decay_t<decltype(arg)>;
        if constexpr (!std::is_same_v<TT, std::shared_ptr<std::monostate>>) {
          SlothInfo::debug("NLSolver::NLSolver setting solver ");
          auto solver = arg.get();
          this->nl_solver_->SetSolver(*solver);
        }
      },
      this->variant_solver_);
  SlothInfo::debug("NLSolver::NLSolver end");
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
  SlothInfo::debug("NLSolver::NLSolver setting operator ");
  this->nl_solver_->SetOperator(ope);

  std::visit(
      [&](auto&& arg) {
        using TT = std::decay_t<decltype(arg)>;
        if constexpr (!std::is_same_v<TT, std::shared_ptr<std::monostate>>) {
          auto solver = arg.get();
          SlothInfo::debug("NLSolver::NLSolver setting solver ");
          this->nl_solver_->SetSolver(*solver);
        }
      },
      this->variant_solver_);
  SlothInfo::debug("NLSolver::NLSolver end");
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

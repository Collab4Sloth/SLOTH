/**
 * @file HPrecondBase.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief  Base class for hypre preconditionners used in linear system  involved in  non linear
 * algorithm
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

class PrecondHypreILU : public SolverBase<mfem::HypreILU, HyprePreconditionerType> {
 public:
  PrecondHypreILU() {}
  std::shared_ptr<mfem::HypreILU> create_solver(HyprePreconditionerType PRECOND,
                                                const Parameters& params) override;

  ~PrecondHypreILU() {}
};

/**
 * @brief Create preconditioner based of the PreconditionerType and a list of Parameters
 *
 * @param PRECOND
 * @param params
 * @return std::shared_ptr<mfem::Solver>
 */
std::shared_ptr<mfem::HypreILU> PrecondHypreILU::create_solver(HyprePreconditionerType PRECOND,
                                                               const Parameters& params) {
  this->solver_description_ = params.get_param_value<std::string>("description");
  SlothInfo::debug(" Create ", this->get_description());

  const int print_level =
      params.get_param_value_or_default<int>("print_level", HYPRE_ILU_DefaultConstant::print_level);
  const HYPRE_Int type = static_cast<HYPRE_Int>(
      params.get_param_value_or_default<int>("type", HYPRE_ILU_DefaultConstant::type));
  const HYPRE_Int n_iter_max = static_cast<HYPRE_Int>(
      params.get_param_value_or_default<int>("iter_max", HYPRE_ILU_DefaultConstant::iter_max));
  const double tol =
      params.get_param_value_or_default<double>("rel_tol", HYPRE_ILU_DefaultConstant::tol);
  const double reordering =
      params.get_param_value_or_default<double>("rel_tol", HYPRE_ILU_DefaultConstant::reorder_type);

  auto pp = std::make_shared<mfem::HypreILU>();
  pp->SetType(type);
  pp->SetMaxIter(n_iter_max);
  pp->SetTol(tol);
  pp->SetLocalReordering(reordering);
  pp->SetPrintLevel(print_level);

  return pp;
}

////////////////////////////////
////////////////////////////////

class PrecondHypreSmoother : public SolverBase<mfem::HypreSmoother, HyprePreconditionerType> {
 public:
  PrecondHypreSmoother() {}
  std::shared_ptr<mfem::HypreSmoother> create_solver(HyprePreconditionerType PRECOND,
                                                     const Parameters& params) override;

  ~PrecondHypreSmoother() {}
};
/**
 * @brief Create preconditioner based of the PreconditionerType and a list of Parameters
 *
 * @param PRECOND
 * @param params
 * @return std::shared_ptr<mfem::Solver>
 */
std::shared_ptr<mfem::HypreSmoother> PrecondHypreSmoother::create_solver(
    HyprePreconditionerType PRECOND, const Parameters& params) {
  this->solver_description_ = params.get_param_value<std::string>("description");
  SlothInfo::debug(" Create ", this->get_description());

  const mfem::HypreSmoother::Type type = static_cast<mfem::HypreSmoother::Type>(
      params.get_param_value_or_default<int>("type", HYPRE_SMOOTHER_DefaultConstant::type));
  const bool positive_diagonal = params.get_param_value_or_default<bool>(
      "positive_diagonal", HYPRE_SMOOTHER_DefaultConstant::positive_diagonal);

  auto pp = std::make_shared<mfem::HypreSmoother>();
  pp->SetType(type);
  pp->SetPositiveDiagonal(positive_diagonal);

  return pp;
}

////////////////////////////////
////////////////////////////////

class PrecondHypreBoomerAMG : public SolverBase<mfem::HypreBoomerAMG, HyprePreconditionerType> {
 public:
  PrecondHypreBoomerAMG() {}
  std::shared_ptr<mfem::HypreBoomerAMG> create_solver(HyprePreconditionerType PRECOND,
                                                      const Parameters& params) override;

  ~PrecondHypreBoomerAMG() {}
};

/**
 * @brief Create preconditioner based of the PreconditionerType and a list of Parameters
 *
 * @param PRECOND
 * @param params
 * @return std::shared_ptr<mfem::Solver>
 */
std::shared_ptr<mfem::HypreBoomerAMG> PrecondHypreBoomerAMG::create_solver(
    HyprePreconditionerType PRECOND, const Parameters& params) {
  this->solver_description_ = params.get_param_value<std::string>("description");
  SlothInfo::debug(" Create ", this->get_description());

  const int print_level = params.get_param_value_or_default<int>(
      "print_level", HYPRE_BOOMER_AMG_DefaultConstant::print_level);
  const HYPRE_Int iter_max = static_cast<HYPRE_Int>(params.get_param_value_or_default<int>(
      "iter_max", HYPRE_BOOMER_AMG_DefaultConstant::iter_max));
  const HYPRE_Real tol = static_cast<HYPRE_Real>(
      params.get_param_value_or_default<double>("iter_max", HYPRE_BOOMER_AMG_DefaultConstant::tol));
  auto pp = std::make_shared<mfem::HypreBoomerAMG>();
  pp->SetMaxIter(iter_max);
  pp->SetPrintLevel(print_level);
  pp->SetTol(tol);

  return pp;
}

////////////////////////////////
////////////////////////////////

class PrecondHypreDiagScale : public SolverBase<mfem::HypreDiagScale, HyprePreconditionerType> {
 public:
  PrecondHypreDiagScale() {}
  std::shared_ptr<mfem::HypreDiagScale> create_solver(HyprePreconditionerType PRECOND,
                                                      const Parameters& params) override;

  ~PrecondHypreDiagScale() {}
};

/**
 * @brief Create preconditioner based of the PreconditionerType and a list of Parameters
 *
 * @param PRECOND
 * @param params
 * @return std::shared_ptr<mfem::HypreDiagScale>
 */
std::shared_ptr<mfem::HypreDiagScale> PrecondHypreDiagScale::create_solver(
    HyprePreconditionerType PRECOND, const Parameters& params) {
  this->solver_description_ = params.get_param_value<std::string>("description");
  SlothInfo::debug(" Create ", this->get_description());

  auto pp = std::make_shared<mfem::HypreDiagScale>();

  return pp;
}

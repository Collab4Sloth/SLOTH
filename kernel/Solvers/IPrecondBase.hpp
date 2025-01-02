/**
 * @file PrecondBase.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief  Base class for preconditionners used in linear system  involved in  non linear algorithm
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

class PrecondDSmoother : public SolverBase<mfem::DSmoother, PreconditionerType> {
 public:
  PrecondDSmoother();
  std::shared_ptr<mfem::DSmoother> create_solver(PreconditionerType PRECOND,
                                                 const Parameters& params) override;

  ~PrecondDSmoother();
};

/**
 * @brief Construct a new PrecondDSmoother::PrecondDSmoother object
 *
 */
PrecondDSmoother::PrecondDSmoother() {}

/**
 * @brief Create preconditioner based of the PreconditionerType and a list of Parameters
 *
 * @param PRECOND
 * @param params
 * @return std::shared_ptr<mfem::Solver>
 */
std::shared_ptr<mfem::DSmoother> PrecondDSmoother::create_solver(PreconditionerType PRECOND,
                                                                 const Parameters& params) {
  // TODO(cci): mettre un check sur le param de base name
  // TODO(cc) : mettre un getinfo pour la doc
  this->solver_description_ = params.get_param_value<std::string>("description");
  SlothInfo::debug(" Create ", this->get_description());

  const int type = params.get_param_value_or_default<int>("type", 1);
  const double scale = params.get_param_value_or_default<double>("scale", 1.);
  const int its = params.get_param_value_or_default<int>("iterations", 1);
  const bool use_abs_diag = params.get_param_value_or_default<int>("use_abs_diag", true);

  auto pp = std::make_shared<mfem::DSmoother>(type, scale, its);
  pp->SetPositiveDiagonal(use_abs_diag);
  return pp;
}

/**
 * @brief Destroy the PrecondDSmoother::PrecondDSmoother object
 *
 */
PrecondDSmoother::~PrecondDSmoother() {}

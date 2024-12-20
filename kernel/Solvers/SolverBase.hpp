
/**
 * @file SolverBase.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Base class for Hypre, Iterative, Direct and Non Linear Solvers
 * @version 0.1
 * @date 2024-07-25
 *
 * Copyright CEA (c) 2024
 *
 */
#include <memory>
#include <string>

#include "Utils/PhaseFieldConstants.hpp"
#include "Utils/PhaseFieldOptions.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

template <class T, class S>
class SolverBase {
 public:
  std::string solver_description_;

  virtual std::shared_ptr<T> create_solver(S SOLVER, const Parameters& params) = 0;

  virtual ~SolverBase() = default;

  std::string get_description() { return this->solver_description_; }
};

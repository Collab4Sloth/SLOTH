/**
 * @file SolverBase.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Base class for Hypre, Iterative, Direct and Non Linear Solvers
 * @version 0.1
 * @date 2025-09-05
 *
 * @copyright 2025 CEA (C)
 *
 * @anchor SolverBase
 *
 *
 * This file is part of SLOTH.
 *
 * SLOTH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SLOTH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#pragma once
#include <memory>
#include <string>

#include "Options/Options.hpp"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

template <class T, class S>
class SolverBase {
 public:
  std::string solver_description_;

  virtual std::shared_ptr<T> create_solver(const Parameters& params) = 0;

  virtual ~SolverBase() = default;

  inline std::string get_description() { return this->solver_description_; }
};

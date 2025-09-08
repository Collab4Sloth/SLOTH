/**
 * @file IPrecondBase.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Base class for preconditionners used in linear system  involved in  non linear algorithm
 * @version 0.1
 * @date 2025-09-05
 * 
 * Copyright CEA (C) 2025
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
#include <memory>
#include <string>

#include "Options/Options.hpp"
#include "Solvers/SolverBase.hpp"
#include "Utils/Utils.hpp"
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

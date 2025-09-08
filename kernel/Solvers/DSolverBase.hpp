/**
 * @file DSolverBase.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Direct solvers
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

class SolverUMFPACK : public SolverBase<mfem::UMFPackSolver, DirectSolverType> {
 public:
  SolverUMFPACK();
  std::shared_ptr<mfem::UMFPackSolver> create_solver(DirectSolverType SOLVER,
                                                     const Parameters& params) override;

  ~SolverUMFPACK();
};

/**
 * @brief Construct a new SolverUMFPACK::SolverUMFPACK object
 *
 */
SolverUMFPACK::SolverUMFPACK() {}

/**
 * @brief Create a direct solver based of the SolverType and a list of Parameters
 *
 * @param SOLVER
 * @param params
 * @return std::shared_ptr<mfem::Solver>
 */
std::shared_ptr<mfem::UMFPackSolver> SolverUMFPACK::create_solver(DirectSolverType SOLVER,
                                                                  const Parameters& params) {
  // TODO(cci): mettre un check sur le param de base name
  // TODO(cc) : mettre un getinfo pour la doc
  this->solver_description_ = params.get_param_value<std::string>("description");
  SlothInfo::debug(" Create ", this->get_description());

  const int print_level =
      params.get_param_value_or_default<int>("print_level", UMFPACK_DefaultConstant::print_level);
  auto ss = std::make_shared<mfem::UMFPackSolver>(MPI_COMM_WORLD);
  ss->SetPrintLevel(print_level);
  return ss;
}

/**
 * @brief Destroy the SolverUMFPACK::SolverUMFPACK object
 *
 */
SolverUMFPACK::~SolverUMFPACK() {}

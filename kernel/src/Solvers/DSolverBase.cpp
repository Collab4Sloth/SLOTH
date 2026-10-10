/**
 * @file DSolverBase.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Direct solvers
 * @version 0.1
 * @date 2025-09-05
 *
 * @copyright CEA (C) 2025
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
#include "Solvers/DSolverBase.hpp"

#include <memory>
#include <string>

#include "Options/Options.hpp"
#include "Parameters/Parameters.hpp"
#include "Solvers/SolverBase.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#ifdef SLOTH_USE_MUMPS
/**
 * @brief Construct a new SolverMUMPS::SolverMUMPS object
 *
 */
SolverMUMPS::SolverMUMPS() {}

/**
 * @brief Create a direct solver based of the SolverType and a list of Parameters
 *
 * @param SOLVER
 * @param params
 * @return std::shared_ptr<mfem::Solver>
 */
std::shared_ptr<mfem::MUMPSSolver> SolverMUMPS::create_solver(const Parameters& params) {
  this->solver_description_ = params.get_param_value<std::string>("description");
  SlothInfo::debug(" Create ", this->get_description());

  int print_level = MUMPS_DefaultConstant::print_level;

  if (verbose_at_least(Verbosity::Verbose)) {
    print_level =
        params.get_param_value_or_default<int>("print_level", MUMPS_DefaultConstant::print_level);
  }

  const int mat_type =
      params.get_param_value_or_default<int>("mat_type", MUMPS_DefaultConstant::mat_type);
  const int reordering_strategy = params.get_param_value_or_default<int>(
      "reordering_strategy", MUMPS_DefaultConstant::reordering_strategy);

  auto ss = std::make_shared<mfem::MUMPSSolver>(MPI_COMM_WORLD);
  ss->SetPrintLevel(print_level);
  ss->SetMatrixSymType(static_cast<mfem::MUMPSSolver::MatType>(mat_type));
  ss->SetReorderingStrategy(
      static_cast<mfem::MUMPSSolver::ReorderingStrategy>(reordering_strategy));
  return ss;
}

/**
 * @brief Destroy the SolverMUMPS::SolverMUMPS object
 *
 */
SolverMUMPS::~SolverMUMPS() {}
#endif
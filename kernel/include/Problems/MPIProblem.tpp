/**
 * @file MPIProblem.tpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Class used to defined a MPI_Problem objet
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

#pragma once
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "Parameters/Parameter.hpp"
#include "PostProcessing/postprocessing.hpp"
#include "Problems/ProblemBase.hpp"
#include "Variables/Variable.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new Unsteady MPI_Problem object
 *
 * @tparam OPE
 * @tparam SOLVER
 * @param oper
 * @param solver
 * @param convergence
 */
template <class VAR, class PST>
MPI_Problem<VAR, PST>::MPI_Problem(VAR& variables, PST& pst)
    : ProblemBase<VAR, PST>("MPI problem", variables, pst) {}

/**
 * @brief Do a time-step by calling Step method of the ODE
 *
 * @tparam OPE
 * @tparam VAR
 * @param unk
 * @param current_time
 * @param current_time_step
 */
template <class VAR, class PST>
void MPI_Problem<VAR, PST>::do_time_step(
    [[maybe_unused]] double& next_time, [[maybe_unused]] const double& current_time,
    [[maybe_unused]] double current_time_step, [[maybe_unused]] const int iter,
    std::vector<std::unique_ptr<mfem::Vector>>& vect_unk,
    [[maybe_unused]] const std::vector<std::vector<std::string>>& unks_info) {
  int rank = mfem::Mpi::WorldRank();

  auto& unk = *(vect_unk[0]);
  unk = static_cast<double>(rank);
  // Store the solution into a temporary mfem::Vector that will be used during updating stage,
  // TODO(cci) : utile? à voir avec le restart de paraview dan mfem et peut etre l'AMR
  this->unknown_.emplace_back(unk);
  // Update time
  next_time = current_time + current_time_step;
}

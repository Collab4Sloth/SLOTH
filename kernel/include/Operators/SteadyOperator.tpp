/**
 * @file SteadyOperator.tpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Steady  Operator
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
#include <vector>

#include "AnalyticalFunctions/AnalyticalFunctions.hpp"
#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Operators/OperatorBase.hpp"
#include "Operators/SteadyReducedOperator.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameters.hpp"
#include "Solvers/LSolver.hpp"
#include "Solvers/NLSolver.hpp"
#include "Spatial/Spatial.hpp"
#include "Utils/Utils.hpp"
#include "Variables/Variables.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new SteadyOperator object
 *
 * @tparam T
 * @tparam DIM
 * @param spatial
 */
template <class T, int DIM>
SteadyOperator<T, DIM>::SteadyOperator(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                                       const std::vector<std::string>& integrators)
    : OperatorBase<T, DIM>(integrators, spatials), steady_reduced_oper(NULL) {
  const Parameters nl_parameters =
      Parameters(Parameter("description", "Newton Algorithm"), Parameter("iterative_mode", true));
  this->overload_nl_solver(NLSolverType::NEWTON, nl_parameters);
}

/**
 * @brief Construct a new SteadyOperator object
 *
 * @tparam T
 * @tparam DIM
 * @param spatial
 * @param source_term_name
 */
template <class T, int DIM>
SteadyOperator<T, DIM>::SteadyOperator(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                                       const std::vector<std::string>& integrators,
                                       std::vector<AnalyticalFunctions<DIM>> source_term_name)
    : OperatorBase<T, DIM>(integrators, spatials, source_term_name), steady_reduced_oper(NULL) {
  const Parameters nl_parameters =
      Parameters(Parameter("description", "Newton Algorithm"), Parameter("iterative_mode", true));
  this->overload_nl_solver(NLSolverType::NEWTON, nl_parameters);
}

/**
 * @brief Construct a new SteadyOperator object
 *
 * @tparam T
 * @tparam DIM
 * @param spatial
 * @param params
 */
template <class T, int DIM>
SteadyOperator<T, DIM>::SteadyOperator(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                                       const std::vector<std::string>& integrators,
                                       const Parameters& params)
    : OperatorBase<T, DIM>(integrators, spatials, params), steady_reduced_oper(NULL) {
  const Parameters nl_parameters =
      Parameters(Parameter("description", "Newton Algorithm"), Parameter("iterative_mode", true));
  this->overload_nl_solver(NLSolverType::NEWTON, nl_parameters);
}

/**
 * @brief Construct a new SteadyOperator object
 *
 * @tparam T
 * @tparam DIM
 * @param spatial
 * @param params
 * @param source_term_name
 */
template <class T, int DIM>
SteadyOperator<T, DIM>::SteadyOperator(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                                       const std::vector<std::string>& integrators,
                                       const Parameters& params,
                                       std::vector<AnalyticalFunctions<DIM>> source_term_name)
    : OperatorBase<T, DIM>(integrators, spatials, params, source_term_name),
      steady_reduced_oper(NULL) {
  const Parameters nl_parameters =
      Parameters(Parameter("description", "Newton Algorithm"), Parameter("iterative_mode", true));
  this->overload_nl_solver(NLSolverType::NEWTON, nl_parameters);
}

/**
 * @brief Initialization stage (call by imeDiscretization<PST, OPE, VAR>::initialize())
 *
 * @tparam T
 * @tparam DIM
 * @param initial_time
 * @param vars
 */
template <class T, int DIM>
void SteadyOperator<T, DIM>::initialize(const double& initial_time, const double time_step,
                                        Variables<T, DIM>& vars,
                                        std::vector<Variables<T, DIM>*> auxvars) {
  Catch_Time_Section("SteadyOperator::initialize");

  OperatorBase<T, DIM>::initialize(initial_time, time_step, vars, auxvars);
}

/**
 * @brief Solve the steady problem
 *
 * @tparam T
 * @tparam DIM
 * @param unk
 * @param current_time
 * @param dt
 */
template <class T, int DIM>
void SteadyOperator<T, DIM>::solve(std::vector<std::unique_ptr<mfem::Vector>>& vect_unk,
                                   double& next_time, const double& current_time, double dt,
                                   const int iter) {
  this->current_time_ = current_time;
  this->current_dt_ = dt;

  // Default iterative_mode to true for Steady Operators
  // Important line if solvers parameters are overloaded.
  this->nl_solver_params_ = this->nl_solver_params_ + Parameter("iterative_mode", true);

  //// Constructing array of offsets
  const size_t unk_size = vect_unk.size();

  //// Constructing BlockVector
  mfem::BlockVector block_unk(this->block_trueOffsets_);
  std::vector<mfem::Vector> u_vect;
  for (size_t i = 0; i < unk_size; i++) {
    auto& unk_i = *(vect_unk[i]);
    mfem::Vector& bb = block_unk.GetBlock(i);
    bb = unk_i;
    u_vect.emplace_back(unk_i);
  }

  this->SetTransientParameters(dt, u_vect);

  /// Apply BCs: check if need to be uncomment
  // for (size_t i = 0; i < unk_size; i++) {
  //   auto &unk_i = *(vect_unk[i]);
  //   this->bcs_[i]->SetBoundaryConditions(unk_i);
  // }

  // Source term
  // const mfem::Array<int> offsets = this->RHS->GetBlockOffsets();
  const int fes_size = this->block_trueOffsets_.Size() - 1;
  mfem::BlockVector source_term(this->block_trueOffsets_);
  source_term = 0.0;
  if (!this->src_func_.empty()) {
    for (int i = 0; i < fes_size; ++i) {
      if (this->src_func_[i] != nullptr) {
        mfem::ParLinearForm* SRC = new mfem::ParLinearForm(this->fes_[i]);
        mfem::Vector& src_i = source_term.GetBlock(i);
        this->get_source_term(i, this->src_func_[i], src_i, SRC);
        delete SRC;
      }
    }
  }
  this->newton_solver_->Mult(source_term, block_unk);

  this->free_memory();

  for (size_t i = 0; i < unk_size; i++) {
    auto& unk_i = *(vect_unk[i]);
    const mfem::Vector& bb = block_unk.GetBlock(i);
    unk_i = bb;
  }
  next_time = current_time + static_cast<double>(iter) * dt;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief Set current dt, unk values - needed to compute action and Jacobian.
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param u_vect unknown vector
 */
template <class T, int DIM>
void SteadyOperator<T, DIM>::SetTransientParameters(const double dt,
                                                    const std::vector<mfem::Vector>& u_vect) {
  Catch_Time_Section("SteadyOperator::SetTransientParameters");

  ////////////////////////////////////////////
  //  Build the RHS of the PDEs
  ////////////////////////////////////////////
  this->build_rhs_nonlinear_form(dt, u_vect);

  ////////////////////////////////////////////
  // Build Newton Linear system
  ////////////////////////////////////////////
  if (this->steady_reduced_oper != nullptr) {
    delete this->steady_reduced_oper;
  }
  this->steady_reduced_oper = new SteadyPhaseFieldReducedOperator(this->RHS, this->ess_tdof_list_);

  ////////////////////////////////////////////
  // Newton Solver
  ////////////////////////////////////////////
  this->SetNewtonAlgorithm(steady_reduced_oper);
}

/**
 * @brief Free memory
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 */
template <class T, int DIM>
void SteadyOperator<T, DIM>::free_memory() {
  delete this->rhs_solver_;
  delete this->RHS;
  this->RHS = nullptr;
  delete this->steady_reduced_oper;
  this->steady_reduced_oper = nullptr;
}

/**
 * @file SteadyOperatorBase.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Steady  Operator base
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

#pragma once

/**
 * @brief SteadyOperatorBase class
 *
 */
template <class T, int DIM>
class SteadyOperatorBase : public OperatorBase<T, DIM> {
 private:
  SteadyPhaseFieldReducedOperator* steady_reduced_oper;

 public:
  SteadyOperatorBase(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                     const std::vector<std::string>& integrators);

  SteadyOperatorBase(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                     const std::vector<std::string>& integrators,
                     std::vector<AnalyticalFunctions<DIM>> source_term_name);

  SteadyOperatorBase(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                     const std::vector<std::string>& integrators, const Parameters& params);

  SteadyOperatorBase(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                     const std::vector<std::string>& integrators, const Parameters& params,
                     std::vector<AnalyticalFunctions<DIM>> source_term_name);

  void Mult([[maybe_unused]] const mfem::Vector& k, [[maybe_unused]] mfem::Vector& y)
      const override {  // Nothing to be done because of manage by steadyreducedoperator
  }
  // mfem::Operator &GetGradient(const mfem::Vector &k) const override;

  // Virtual methods

  // void initialize(const double &initial_time, Variables<T, DIM> &vars) override;
  void initialize(const double& initial_time, Variables<T, DIM>& vars,
                  std::vector<Variables<T, DIM>*> auxvars) override;
  // Pure virtual methods
  void SetTransientParameters(const std::vector<mfem::Vector>& u_vet) override;
  void solve(std::vector<std::unique_ptr<mfem::Vector>>& vect_unk, double& next_time,
             const double& current_time, double dt, const int iter) override;
  SlothNLFormIntegrator<Variables<T, DIM>>* set_nlfi_ptr(
      const std::string nlfi, const std::vector<mfem::Vector>& u) override = 0;
  void get_parameters() override = 0;
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Construct a new SteadyOperatorBase object
 *
 * @tparam T
 * @tparam DIM
 * @param spatial
 */
template <class T, int DIM>
SteadyOperatorBase<T, DIM>::SteadyOperatorBase(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                                               const std::vector<std::string>& integrators)
    : OperatorBase<T, DIM>(integrators, spatials), steady_reduced_oper(NULL) {
  const Parameters nl_parameters =
      Parameters(Parameter("description", "Newton Algorithm"), Parameter("iterative_mode", true));
  this->overload_nl_solver(NLSolverType::NEWTON, nl_parameters);
}

/**
 * @brief Construct a new SteadyOperatorBase object
 *
 * @tparam T
 * @tparam DIM
 * @param spatial
 * @param source_term_name
 */
template <class T, int DIM>
SteadyOperatorBase<T, DIM>::SteadyOperatorBase(
    std::vector<SpatialDiscretization<T, DIM>*> spatials,
    const std::vector<std::string>& integrators,
    std::vector<AnalyticalFunctions<DIM>> source_term_name)
    : OperatorBase<T, DIM>(integrators, spatials, source_term_name), steady_reduced_oper(NULL) {
  const Parameters nl_parameters =
      Parameters(Parameter("description", "Newton Algorithm"), Parameter("iterative_mode", true));
  this->overload_nl_solver(NLSolverType::NEWTON, nl_parameters);
}

/**
 * @brief Construct a new SteadyOperatorBase object
 *
 * @tparam T
 * @tparam DIM
 * @param spatial
 * @param params
 */
template <class T, int DIM>
SteadyOperatorBase<T, DIM>::SteadyOperatorBase(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                                               const std::vector<std::string>& integrators,
                                               const Parameters& params)
    : OperatorBase<T, DIM>(integrators, spatials, params), steady_reduced_oper(NULL) {
  const Parameters nl_parameters =
      Parameters(Parameter("description", "Newton Algorithm"), Parameter("iterative_mode", true));
  this->overload_nl_solver(NLSolverType::NEWTON, nl_parameters);
}

/**
 * @brief Construct a new SteadyOperatorBase object
 *
 * @tparam T
 * @tparam DIM
 * @param spatial
 * @param params
 * @param source_term_name
 */
template <class T, int DIM>
SteadyOperatorBase<T, DIM>::SteadyOperatorBase(
    std::vector<SpatialDiscretization<T, DIM>*> spatials,
    const std::vector<std::string>& integrators, const Parameters& params,
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
void SteadyOperatorBase<T, DIM>::initialize(const double& initial_time, Variables<T, DIM>& vars,
                                            std::vector<Variables<T, DIM>*> auxvars) {
  Catch_Time_Section("SteadyOperatorBase::initialize");

  OperatorBase<T, DIM>::initialize(initial_time, vars, auxvars);
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
void SteadyOperatorBase<T, DIM>::solve(std::vector<std::unique_ptr<mfem::Vector>>& vect_unk,
                                       double& next_time, const double& current_time, double dt,
                                       const int iter) {
  this->current_time_ = current_time;
  this->current_dt_ = dt;

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

  this->SetTransientParameters(u_vect);

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
        mfem::ParLinearForm* RHS = new mfem::ParLinearForm(this->fes_[i]);
        mfem::Vector& src_i = source_term.GetBlock(i);
        this->get_source_term(i, this->src_func_[i], src_i, RHS);
        delete RHS;
      }
    }
  }
  this->newton_solver_->Mult(source_term, block_unk);
  delete this->rhs_solver_;

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
 *solution_coef
 * @param dt time-step
 * @param u unknown vector
 */
template <class T, int DIM>
void SteadyOperatorBase<T, DIM>::SetTransientParameters(const std::vector<mfem::Vector>& u_vect) {
  Catch_Time_Section("SteadyOperatorBase::SetTransientParameters");

  ////////////////////////////////////////////
  //  Build the RHS of the PDEs
  ////////////////////////////////////////////
  this->build_rhs_nonlinear_form(u_vect);

  ////////////////////////////////////////////
  // Build Newton Linear system
  ////////////////////////////////////////////
  if (steady_reduced_oper != nullptr) {
    delete steady_reduced_oper;
  }
  steady_reduced_oper = new SteadyPhaseFieldReducedOperator(this->RHS, this->ess_tdof_list_);

  ////////////////////////////////////////////
  // Newton Solver
  ////////////////////////////////////////////
  this->SetNewtonAlgorithm(steady_reduced_oper);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

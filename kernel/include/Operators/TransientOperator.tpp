/**
 * @file TransientOperator.tpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief  Transient Operator
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
#include <algorithm>
#include <cmath>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "AnalyticalFunctions/AnalyticalFunctions.hpp"
#include "Coefficients/MfemCoefficient.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Operators/OperatorBase.hpp"
#include "Operators/ReducedOperator.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Solvers/LSolver.hpp"
#include "Spatial/Spatial.hpp"
#include "Time/TimeStepFailure.hpp"
#include "Utils/Utils.hpp"
#include "Variables/Variables.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Construct a new TransientOperator object
 *
 * @tparam T
 * @tparam DIM
 * @param spatial
 * @param ode
 */
template <class T, int DIM>
TransientOperator<T, DIM>::TransientOperator(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                                             const std::vector<std::string>& rhs_integrators,
                                             TimeScheme::value ode,
                                             const std::string lhs_integrator)
    : OperatorBase<T, DIM>(rhs_integrators, spatials),
      mfem::TimeDependentOperator(this->compute_total_height(spatials),
                                  this->compute_total_width(spatials), 0.0),
      lhs_integrator_(lhs_integrator),
      M(NULL),
      LHS(NULL),
      reduced_oper(NULL) {
  this->set_ODE_solver(ode);
  this->set_default_mass_solver();
}

/**
 * @brief Construct a new TransientOperator object
 *
 * @tparam T
 * @tparam DIM
 * @param spatial
 * @param params
 * @param ode
 */
template <class T, int DIM>
TransientOperator<T, DIM>::TransientOperator(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                                             const std::vector<std::string>& rhs_integrators,
                                             const Parameters& params, TimeScheme::value ode,
                                             const std::string lhs_integrator)
    : OperatorBase<T, DIM>(rhs_integrators, spatials, params),
      mfem::TimeDependentOperator(this->compute_total_height(spatials),
                                  this->compute_total_width(spatials), 0.0),
      lhs_integrator_(lhs_integrator),
      M(NULL),
      LHS(NULL),
      reduced_oper(NULL) {
  this->set_ODE_solver(ode);
  this->set_default_mass_solver();
}

/**
 * @brief Destroy the TransientOperator object
 *        Solvers are released before the matrices they refer to.
 *
 * @tparam T
 * @tparam DIM
 */
template <class T, int DIM>
TransientOperator<T, DIM>::~TransientOperator() {
  this->free_mass_matrices();
  delete this->M;
  this->M = nullptr;
  delete this->reduced_oper;
  this->reduced_oper = nullptr;
  delete this->LHS;
  this->LHS = nullptr;
}

/**
 * @brief  Set the ODE time marching
 *
 * @tparam T
 * @tparam DIM
 * @param ode_solver
 */
template <class T, int DIM>
void TransientOperator<T, DIM>::set_ODE_solver(const TimeScheme::value& ode_solver) {
  // TODO(cci) faire un template?
  switch (ode_solver) {
    case TimeScheme::EulerExplicit: {
      this->ode_solver_ = std::make_shared<mfem::ForwardEulerSolver>();
      this->is_explicit_ = true;
      break;
    }
    case TimeScheme::EulerImplicit: {
      this->ode_solver_ = std::make_shared<mfem::BackwardEulerSolver>();
      this->is_explicit_ = false;
      break;
    }
    case TimeScheme::RungeKutta4: {
      // explicit forth-order Runge-Kutta
      this->ode_solver_ = std::make_shared<mfem::RK4Solver>();
      this->is_explicit_ = true;
      break;
    }
    case TimeScheme::SDIRK23: {
      //  two-step third-order Singly-diagonal implicit Runge-Kutta scheme (SDIRK23)
      this->ode_solver_ = std::make_shared<mfem::SDIRK23Solver>();
      this->is_explicit_ = false;
      break;
    }
    case TimeScheme::SDIRK33: {
      //  three-step third-order Singly-diagonal implicit Runge-Kutta scheme (SDIRK23)
      this->ode_solver_ = std::make_shared<mfem::SDIRK33Solver>();
      this->is_explicit_ = false;
      break;
    }
    case TimeScheme::BDF2: {
      // variable-step second-order backward differentiation formula (first step: implicit Euler)
      this->bdf2_solver_ = std::make_shared<SlothBDF2Solver<T, DIM>>();
      this->ode_solver_ = this->bdf2_solver_;
      this->is_explicit_ = false;
      break;
    }
    default:
      mfem::mfem_error(
          "TimeDiscretization::set_ODE_solver: EulerImplicit, EulerExplicit, RungeKutta4, SDIRK23, "
          "SDIRK33, BDF2 are available");
      break;
  }
}

/**
 * @brief  Initialization stage (call by imeDiscretization<PST, OPE, VAR>::initialize())
 *
 * @tparam T
 * @tparam DIM
 * @param initial_time
 * @param vars
 */
template <class T, int DIM>
void TransientOperator<T, DIM>::initialize(const double& initial_time, const double time_step,
                                           Variables<T, DIM>& vars,
                                           std::vector<Variables<T, DIM>*> auxvars) {
  Catch_Time_Section("TransientOperator::initialize");

  this->SetTime(initial_time);

  this->vars_ = &vars;

  OperatorBase<T, DIM>::initialize(initial_time, time_step, vars, auxvars);

  if (this->bdf2_solver_ != nullptr) {
    this->bdf2_solver_->set_variables(vars, auxvars);
  }
  this->ode_solver_->Init(*this);

  // Get coefficients for time-derivative in case of explicit solvers
  if (this->is_explicit_) {
    this->get_explicit_time_coefficients();
  }
}

/**
 * @brief Solve the time-dependent problem
 *
 * @tparam T
 * @tparam DIM
 * @param unk
 * @param current_time
 * @param current_time_step
 */
template <class T, int DIM>
void TransientOperator<T, DIM>::solve(std::vector<std::unique_ptr<mfem::Vector>>& vect_unk,
                                      double& next_time, const double& current_time,
                                      double current_time_step, [[maybe_unused]] const int iter) {
  if (this->enable_amr_ && this->mesh_has_changed()) {
    // Resynchronize height/width/z/block_trueOffsets_
    this->UpdateAfterMeshChange();
    // Resynchronize the OTHER height/width copy, inherited via mfem::TimeDependentOperator,
    // which UpdateAfterMeshChange() does not touch.
    this->mfem::TimeDependentOperator::height = this->OperatorBase<T, DIM>::height;
    this->mfem::TimeDependentOperator::width = this->OperatorBase<T, DIM>::width;
    // Resize k (the ODE solver's internal work vector)
    this->ode_solver_->Init(*this);
    // The mesh has changed: the mass matrices must be rebuilt
    this->free_mass_matrices();
    // Split systems: lumped masses and last projected state refer to the old mesh
    this->lumped_.clear();
    this->last_projection_signature_.Destroy();
    this->free_constraint_solvers();
    this->record_fes_sequence();
  }
  //// Constructing array of offsets
  const size_t unk_size = vect_unk.size();
  this->solve_time_step_ = current_time_step;
  //// Constructing BlockVector
  mfem::BlockVector block_unk(this->block_trueOffsets_);
  for (size_t i = 0; i < unk_size; i++) {
    auto& unk_i = *(vect_unk[i]);
    mfem::Vector& bb = block_unk.GetBlock(static_cast<int>(i));
    bb = unk_i;
  }
  //// Other schemes: lagged auxiliary variables of a partitioned coupling reported once (BDF2
  //// extrapolates them, see SlothBDF2Solver)
  this->report_lagged_auxiliary_variables();

  //// Call ODE solver
  this->current_time_ = current_time;
  this->current_dt_ = current_time_step;
  this->newton_iterations_ = -1;
  // Throws TimeStepFailure if a nonlinear solve fails (see ImplicitSolve): vect_unk is then left
  // unchanged
  this->ode_solver_->Step(block_unk, next_time, current_time_step);

  // The last stage of a Runge-Kutta scheme is not necessarily at the end of the step (e.g.
  // SDIRK23): coefficients at the end of the step for the post-processing
  this->SetTime(next_time);

  // Split (DAE) systems with explicit schemes: the algebraic variables are not advanced by the
  // ODE solver (d/dt = 0). They are made consistent with the differential variables at the end
  // of the step (post-processing, next step)
  if (this->split_explicit_) {
    std::vector<mfem::Vector> views;
    views.reserve(unk_size);
    for (size_t i = 0; i < unk_size; i++) {
      mfem::Vector& bi = block_unk.GetBlock(static_cast<int>(i));
      views.emplace_back(bi.GetData(), bi.Size());
    }
    this->build_rhs_nonlinear_form(current_time_step, views);
    this->project_algebraic_variables(block_unk);
  }

  //// Updating vect_unk
  for (size_t i = 0; i < unk_size; i++) {
    auto& unk_i = *(vect_unk[i]);
    const mfem::Vector& bb = block_unk.GetBlock(static_cast<int>(i));
    unk_i = bb;
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief Partitioned coupling, schemes other than BDF2: report once the lagged auxiliary variables
 *
 * @details An auxiliary variable belonging to a problem not yet solved in the current step is
 *          seen at t^n (lag O(dt), first order) and its previous value is that of the previous
 *          step. BDF2 extrapolates these variables instead (SlothBDF2Solver).
 *
 * @tparam T
 * @tparam DIM
 */
template <class T, int DIM>
void TransientOperator<T, DIM>::report_lagged_auxiliary_variables() {
  if (this->bdf2_solver_ != nullptr || this->lagged_aux_info_printed_ || this->vars_ == nullptr) {
    return;
  }
  if (this->vars_->get_variables_number() == 0 ||
      (*this->vars_)[0].get_number_of_updates() == 0) {
    return;  // first step: fixed fields cannot be distinguished
  }
  this->lagged_aux_info_printed_ = true;
  if (mfem::Mpi::WorldRank() != 0) return;
  for (auto* auxvar_vec : this->auxvariables_) {
    if (auxvar_vec == nullptr) continue;
    for (size_t i = 0; i < auxvar_vec->get_variables_number(); i++) {
      auto& auxvar = (*auxvar_vec)[i];
      if (is_lagged_auxiliary_variable(auxvar, *this->vars_)) {
        SlothInfo::verbose(" TransientOperator: auxiliary variable '", auxvar.getVariableName(),
                           "' not yet advanced in the step (partitioned coupling): seen at the "
                           "beginning of the step (lag O(dt)), its previous value is that of the "
                           "previous step");
      }
    }
  }
}

/**
 * @brief Release the mass solvers, then the mass matrices they refer to
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 */
template <class T, int DIM>
void TransientOperator<T, DIM>::free_mass_matrices() {
  // Order matters: solvers first, matrices afterwards
  this->M_solver_.clear();
  this->mass_matrix_solver_.reset();
  this->mass_solvers_.clear();
  for (auto* mat : this->Mmat_) {
    delete mat;
  }
  this->Mmat_.clear();
}

/**
 * @brief build the mass matrix (one lumped mass matrix and one solver per block)
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 */
template <class T, int DIM>
void TransientOperator<T, DIM>::build_mass_matrix(const std::vector<mfem::Vector>& u_vect) {
  this->free_mass_matrices();

  std::vector<mfem::ParGridFunction> vauxn;
  for (const auto& auxvar_vec : this->auxvariables_) {
    for (auto auxvars : auxvar_vec->getVariables()) {
      auto fes = auxvars.get_fespace();
      mfem::ParGridFunction auxn(fes);
      auto auxvar_n = auxvars.get_second_to_last();
      auxn.SetFromTrueDofs(auxvar_n);
      vauxn.emplace_back(auxn);
    }
  }

  const size_t nblocks = u_vect.size();
  this->Mmat_.reserve(nblocks);
  this->mass_solvers_.reserve(nblocks);
  this->M_solver_.reserve(nblocks);

  for (size_t i = 0; i < nblocks; i++) {
    const int ib = static_cast<int>(i);
    if (M != nullptr) {
      delete M;
    }
    M = new mfem::ParBilinearForm(this->fes_[ib]);

    // Build Mfem Coefficient from Sloth Coefficient
    // Sloth coefficients are either scalar or explicit (variable + {auxiliary variables})
    mfem::ParGridFunction un(this->fes_[ib]);
    un.SetFromTrueDofs(u_vect[i]);

    auto coef_a_exp = MfemCoefficient(this->explicit_time_coefficients_[0], {un}, {un}, vauxn);
    auto coef_b_exp = MfemCoefficient(this->explicit_time_coefficients_[1], {un}, {un}, vauxn);
    mfem::ProductCoefficient mass_coefficient(coef_a_exp, coef_b_exp);

    M->AddDomainIntegrator(new mfem::LumpedIntegrator(new mfem::MassIntegrator(mass_coefficient)));
    M->Assemble(0);
    M->Finalize(0);

    this->Mmat_.emplace_back(M->ParallelAssemble());
    std::unique_ptr<mfem::HypreParMatrix> Me(
        this->Mmat_[i]->EliminateRowsCols(this->ess_tdof_list_[i]));

    this->mass_solvers_.emplace_back(std::make_shared<LSolver>(
        this->mass_solver_, this->mass_solver_params_, this->mass_precond_,
        this->mass_precond_params_, *this->Mmat_[i]));
    this->M_solver_.emplace_back(this->mass_solvers_.back()->get_solver());
  }

  // The form M holds integrators referring to local coefficients: release it now
  delete M;
  M = nullptr;
}

/**
 * @brief Build the nonlinear form associated with LHS of the PDEs
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param dt
 * @param u
 */
template <class T, int DIM>
void TransientOperator<T, DIM>::build_lhs_nonlinear_form(const double dt,
                                                         const std::vector<mfem::Vector>& u_vect) {
  if (LHS != nullptr) {
    delete LHS;
  }
  LHS = new mfem::ParBlockNonlinearForm(this->fes_);
  auto integrator_ptr = this->set_lhs_nlfi_ptr(dt, u_vect);
  this->LHS->AddDomainIntegrator(integrator_ptr);
}

/**
 * @brief  Set the LHS NonLinearFormIntegrator
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param dt
 * @param u_vect
 */
template <class T, int DIM>
SlothNLFormIntegrator<Variables<T, DIM>>* TransientOperator<T, DIM>::set_lhs_nlfi_ptr(
    const double dt, const std::vector<mfem::Vector>& u) {
  Catch_Time_Section("TransientOperator::set_lhs_nlfi_ptr");

  std::vector<mfem::ParGridFunction> vun;
  std::vector<mfem::ParGridFunction> vauxn;
  for (size_t i = 0; i < u.size(); i++) {
    mfem::ParGridFunction un(this->fes_[static_cast<int>(i)]);
    un.SetFromTrueDofs(u[i]);
    vun.emplace_back(un);
  }
  for (const auto& auxvar_vec : this->auxvariables_) {
    for (auto auxvars : auxvar_vec->getVariables()) {
      auto fes = auxvars.get_fespace();
      mfem::ParGridFunction auxn(fes);
      auto auxvar_n = auxvars.get_second_to_last();
      auxn.SetFromTrueDofs(auxvar_n);
      vauxn.emplace_back(auxn);
    }
  }
  const Parameters& all_params = this->params_ - this->default_p_;

  auto lhs = this->get_lhs_integrator(dt, this->lhs_integrator_, vun, vauxn, all_params);
  lhs->set_physical_time_step(this->current_dt_ > 0. ? this->current_dt_ : dt);
  lhs->init();
  return lhs;
}

/**
 * @brief Set current dt, unk values - needed to compute action and Jacobian.
 *solution_coef
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param dt time-step
 * @param u unknown vector
 */
template <class T, int DIM>
void TransientOperator<T, DIM>::SetTransientParameters(const double dt,
                                                       const std::vector<mfem::Vector>& u_vect) {
  Catch_Time_Section("TransientOperator::SetTransientParameters");

  ////////////////////////////////////////////
  // Build the LHS of the PDEs
  ////////////////////////////////////////////
  this->build_lhs_nonlinear_form(dt, u_vect);

  ////////////////////////////////////////////
  //  Build the RHS of the PDEs
  ////////////////////////////////////////////
  this->build_rhs_nonlinear_form(dt, u_vect);
  ////////////////////////////////////////////
  // Build Newton Linear system
  ////////////////////////////////////////////
  if (this->reduced_oper != nullptr) {
    delete this->reduced_oper;
  }
  this->reduced_oper = new PhaseFieldReducedOperator(this->LHS, this->RHS, this->ess_tdof_list_);
  ////////////////////////////////////////////
  // Newton Solver
  ////////////////////////////////////////////
  this->SetNewtonAlgorithm(this->reduced_oper);
}
/**
 * @brief Compute the mass matrix and the non linear form with the solution at the previous time
 * step
 *
 * @details For split (DAE) systems ("Split" in the name of the LHS integrator), the mass is
 *          taken from the LHS integrator itself (lumped Jacobian blocks), so that the explicit
 *          and implicit paths use exactly the same time-derivative terms and coefficients.
 *
 * @tparam T
 * @tparam DIM
 */
template <class T, int DIM>
void TransientOperator<T, DIM>::SetExplicitTransientParameters(
    const double dt, const std::vector<mfem::Vector>& un_vect) {
  Catch_Time_Section("TransientOperator::SetExplicitTransientParameters");

  ////////////////////////////////////////////
  // Structure of the system (once)
  ////////////////////////////////////////////
  if (!this->split_checked_) {
    this->setup_split_structure(dt, un_vect);
  }

  ////////////////////////////////////////////
  // Mass matrix
  ////////////////////////////////////////////
  if (this->split_explicit_) {
    // Mass taken from the LHS integrator (same coefficients as the implicit path).
    // Built once if every ExplicitTime coefficient is scalar, otherwise at each call.
    if (!this->split_mass_constant_ || this->lumped_.size() != this->evol_.size()) {
      this->build_lhs_nonlinear_form(dt, un_vect);
      this->compute_split_lumped_masses();
    }
  } else if (!this->constant_mass_matrix_ || this->M_solver_.size() != un_vect.size()) {
    // built once if constant, otherwise at each call
    this->build_mass_matrix(un_vect);
  }

  ////////////////////////////////////////////
  // PhaseField non linear form
  ////////////////////////////////////////////
  this->build_rhs_nonlinear_form(dt, un_vect);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief Apply a SLOTH linear solver (variant) to a vector
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param solver SLOTH solver (variant)
 * @param in input vector
 * @param out output vector
 */
template <class T, int DIM>
void TransientOperator<T, DIM>::apply_solver(const VSharedMFEMSolver& solver,
                                             const mfem::Vector& in, mfem::Vector& out) {
  std::visit(
      [&](auto&& arg) {
        using TT = std::decay_t<decltype(arg)>;
        if constexpr (!std::is_same_v<TT, std::shared_ptr<std::monostate>>) {
          if (in.Size() != arg->NumCols()) {
            std::string msg = "Size mismatch: in.Size() [" + std::to_string(in.Size()) +
                              "] != arg->NumCols() [" + std::to_string(arg->NumCols()) + "]";
            MFEM_ABORT(msg.c_str());
          }
          arg->Mult(in, out);
        } else {
          MFEM_ABORT("TransientOperator::apply_solver: no solver available.");
        }
      },
      solver);
}

/**
 * @brief Analyse the structure of a split (DAE) system solved with an explicit scheme.
 *
 * @details Active if the scheme is explicit and the name of the LHS integrator contains
 *          "Split". The LHS integrator gives the (row, variable) pairs such that row 'row'
 *          carries d(variable)/dt (time_derivative_map()). The other rows are algebraic
 *          constraints, solved for the variables without time derivative (paired in order).
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param dt time step
 * @param u_vect unknown vectors
 */
template <class T, int DIM>
void TransientOperator<T, DIM>::setup_split_structure(const double dt,
                                                      const std::vector<mfem::Vector>& u_vect) {
  this->split_checked_ = true;
  this->split_explicit_ =
      this->is_explicit_ && this->lhs_integrator_.find("Split") != std::string::npos;
  if (!this->split_explicit_) return;

  std::unique_ptr<SlothNLFormIntegrator<Variables<T, DIM>>> lhs(this->set_lhs_nlfi_ptr(dt, u_vect));
  const auto map = lhs->time_derivative_map();
  MFEM_VERIFY(!map.empty(), "Split LHS integrator '" + this->lhs_integrator_ +
                                "' must define time_derivative_map() for explicit schemes.");

  const int n = this->fes_.Size();
  std::vector<bool> row_evol(static_cast<size_t>(n), false);
  std::vector<bool> var_evol(static_cast<size_t>(n), false);
  this->evol_.clear();
  this->alg_.clear();
  for (const auto& [r, c] : map) {
    MFEM_VERIFY(r >= 0 && r < n && c >= 0 && c < n,
                "Split system: time_derivative_map() refers to an unknown block.");
    this->evol_.emplace_back(r, c);
    row_evol[static_cast<size_t>(r)] = true;
    var_evol[static_cast<size_t>(c)] = true;
  }
  std::vector<int> alg_rows;
  std::vector<int> alg_vars;
  for (int i = 0; i < n; ++i) {
    if (!row_evol[static_cast<size_t>(i)]) alg_rows.push_back(i);
    if (!var_evol[static_cast<size_t>(i)]) alg_vars.push_back(i);
  }
  MFEM_VERIFY(alg_rows.size() == alg_vars.size(),
              "Split system: number of constraint rows != number of algebraic variables.");
  for (size_t k = 0; k < alg_rows.size(); ++k) {
    this->alg_.emplace_back(alg_rows[k], alg_vars[k]);
  }

  for (const auto& [r, c] : this->evol_) {
    MFEM_VERIFY(this->fes_[r]->GetTrueVSize() == this->fes_[c]->GetTrueVSize(),
                "Split system: row and variable spaces must have the same size.");
  }
  for (const auto& [r, c] : this->alg_) {
    MFEM_VERIFY(this->fes_[r]->GetTrueVSize() == this->fes_[c]->GetTrueVSize(),
                "Split system: row and variable spaces must have the same size.");
    MFEM_VERIFY(this->ess_tdof_list_[static_cast<size_t>(c)].Size() == 0,
                "Split system: Dirichlet conditions on an algebraic variable are not supported.");
  }

  // Modified Newton of the projection: minimal residual reduction per correction before the
  // constraint Jacobian is rebuilt (optional parameter)
  this->constraint_min_reduction_ = this->params_.template get_param_value_or_default<double>(
      "constraint_jacobian_min_reduction", 0.1);
  MFEM_VERIFY(this->constraint_min_reduction_ > 0. && this->constraint_min_reduction_ < 1.,
              "constraint_jacobian_min_reduction must be in ]0, 1[. Please check your data.");

  // Lumped masses built once if every ExplicitTime coefficient (read by the LHS integrator for
  // each block) is scalar or absent; otherwise rebuilt at each stage
  this->split_mass_constant_ = true;
  for (int b = 0; b < n; ++b) {
    auto coef = this->get_coefficient(b, GlossaryType::ExplicitTime, 0);
    if (coef.has_value() && !(*coef).is_scalar()) {
      this->split_mass_constant_ = false;
    }
  }
}

/**
 * @brief Lumped masses of a split system, taken from the LHS integrator.
 *
 * @details The LHS Jacobian is evaluated at k = 0 (time derivative): the coefficients are then
 *          evaluated with the state of the stage and the terms in dt dA/du vanish. Row-sum
 *          lumping of the block (row, variable) of each evol_ pair.
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 */
template <class T, int DIM>
void TransientOperator<T, DIM>::compute_split_lumped_masses() {
  mfem::BlockVector k0(this->block_trueOffsets_);
  k0 = 0.;
  auto& JL = this->LHS->GetGradient(k0);
  this->lumped_.resize(this->evol_.size());
  for (size_t p = 0; p < this->evol_.size(); ++p) {
    const auto [r, c] = this->evol_[p];
    mfem::Operator& B = JL.GetBlock(r, c);
    mfem::Vector ones(B.Width());
    ones = 1.;
    this->lumped_[p].SetSize(B.Height());
    B.Mult(ones, this->lumped_[p]);  // row-sum lumping
  }
}

/**
 * @brief Make the algebraic variables of a split system consistent with the differential ones.
 *
 * @details Modified Newton iterations on the constraint rows with respect to the algebraic
 *          variables. The constraint Jacobian blocks (and their solvers) are kept between
 *          projections and refreshed at the current state only if they are missing or if the
 *          last correction did not reduce the residual by at least a factor 10. Correctness
 *          relies on the residual test only: a stale Jacobian can only cost iterations.
 *          For constraints affine in the algebraic variables with a constant Jacobian (e.g.
 *          Cahn-Hilliard, split Allen-Cahn), one correction without any Jacobian assembly.
 *          The linear systems use the solver settings of the mass matrix.
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param ub unknown block vector (algebraic blocks are updated in place)
 */
template <class T, int DIM>
void TransientOperator<T, DIM>::project_algebraic_variables(mfem::BlockVector& ub) {
  Catch_Time_Section("TransientOperator::project_algebraic_variables");
  if (this->alg_.empty()) return;

  MPI_Comm comm = this->fes_[0]->GetComm();
  mfem::BlockVector res(this->block_trueOffsets_);
  double norm0 = -1.;
  double norm_prev = -1.;
  const int max_iter = 20;
  const double min_reduction = 0.1;  // refresh the Jacobian if norm > 0.1 * previous norm
  for (int it = 0; it < max_iter; ++it) {
    this->RHS->Mult(ub, res);
    double norm2 = 0.;
    for (const auto& [r, c] : this->alg_) {
      norm2 += mfem::InnerProduct(comm, res.GetBlock(r), res.GetBlock(r));
    }
    const double norm = std::sqrt(norm2);
    if (norm0 < 0.) norm0 = norm;
    if (norm <= 1.e-14 || norm <= 1.e-10 * norm0) {
      // Remember the consistent state (skipped if passed again unchanged, see Mult)
      this->last_projection_signature_ = this->projection_signature(ub);
      return;
    }

    // Jacobian blocks: missing, or too slow convergence with the kept ones -> rebuild here
    const bool missing = this->constraint_solvers_.size() != this->alg_.size();
    const bool slow = norm_prev > 0. && norm > min_reduction * norm_prev;
    if (missing || slow) {
      this->build_constraint_solvers(ub);
    }
    norm_prev = norm;

    for (size_t k = 0; k < this->alg_.size(); ++k) {
      const auto [r, c] = this->alg_[k];
      mfem::Vector d(ub.GetBlock(c).Size());
      d = 0.;
      apply_solver(this->constraint_solvers_[k]->get_solver(), res.GetBlock(r), d);
      ub.GetBlock(c) -= d;
    }
  }
  // Same decision on every rank (norms reduced over all the ranks)
  throw TimeStepFailure("Split system: projection of the algebraic variables did not converge.");
}

/**
 * @brief Assemble the constraint Jacobian blocks at the state ub, copy them and build their
 *        solvers (solver settings of the mass matrix).
 *
 * @details Deep copies: the gradient belongs to the RHS form, which is rebuilt at each stage.
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param ub unknown block vector
 */
template <class T, int DIM>
void TransientOperator<T, DIM>::build_constraint_solvers(const mfem::BlockVector& ub) {
  Catch_Time_Section("TransientOperator::build_constraint_solvers");
  this->free_constraint_solvers();
  auto& J = this->RHS->GetGradient(ub);
  for (const auto& [r, c] : this->alg_) {
    auto* Jrc = dynamic_cast<mfem::HypreParMatrix*>(&J.GetBlock(r, c));
    MFEM_VERIFY(Jrc != nullptr, "Split system: constraint Jacobian block is not a HypreParMatrix.");
    this->constraint_mats_.emplace_back(std::make_shared<mfem::HypreParMatrix>(*Jrc));
    this->constraint_solvers_.emplace_back(std::make_shared<LSolver>(
        this->mass_solver_, this->mass_solver_params_, this->mass_precond_,
        this->mass_precond_params_, *this->constraint_mats_.back()));
  }
}

/**
 * @brief Release the constraint solvers, then the matrices they refer to.
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 */
template <class T, int DIM>
void TransientOperator<T, DIM>::free_constraint_solvers() {
  // Order matters: solvers first, matrices afterwards
  this->constraint_solvers_.clear();
  this->constraint_mats_.clear();
}

/**
 * @brief Signature of a state for the projection of the algebraic variables.
 *
 * @details Everything the constraint rows may depend on: the unknowns, the auxiliary variables
 *          (as read by the integrators), the time (time-dependent coefficients) and the time
 *          step. Two identical signatures give the same projection.
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param u unknown vector
 * @return mfem::Vector signature
 */
template <class T, int DIM>
mfem::Vector TransientOperator<T, DIM>::projection_signature(const mfem::Vector& u) const {
  std::vector<mfem::Vector> aux;
  int size = u.Size() + 2;
  for (const auto& auxvar_vec : this->auxvariables_) {
    for (auto auxvar : auxvar_vec->getVariables()) {
      aux.emplace_back(auxvar.get_second_to_last());
      size += aux.back().Size();
    }
  }
  mfem::Vector sig(size);
  int k = 0;
  for (int i = 0; i < u.Size(); ++i) sig[k++] = u[i];
  for (const auto& a : aux) {
    for (int i = 0; i < a.Size(); ++i) sig[k++] = a[i];
  }
  sig[k++] = this->mfem::TimeDependentOperator::GetTime();
  sig[k++] = this->solve_time_step_;
  return sig;
}

/**
 * @brief True if the state u has already been projected (identical signature, bitwise), on all
 *        MPI ranks (collective decision: the projection involves collective operations).
 *
 * @details Typical case: the first stage of an explicit Runge-Kutta scheme is evaluated at the
 *          state of the end of the previous step, already projected in solve().
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param u unknown vector
 */
template <class T, int DIM>
bool TransientOperator<T, DIM>::is_already_projected(const mfem::Vector& u) const {
  int local_same = 0;
  if (this->last_projection_signature_.Size() > 0) {
    const mfem::Vector sig = this->projection_signature(u);
    local_same = (sig.Size() == this->last_projection_signature_.Size()) ? 1 : 0;
    for (int k = 0; local_same && k < sig.Size(); ++k) {
      if (sig[k] != this->last_projection_signature_[k]) local_same = 0;
    }
  }
  int global_same = 0;
  MPI_Allreduce(&local_same, &global_same, 1, MPI_INT, MPI_LAND, this->fes_[0]->GetComm());
  return global_same == 1;
}

/**
 * @brief Compute the right-hand side of the ODE system.
 *
 * @details Standard systems: M du_i/dt = -R_i for each block i.
 *          Split (DAE) systems with explicit schemes:
 *            1) the algebraic variables are made consistent with the differential ones,
 *            2) lumped(row, var) d(var)/dt = -R_row for each (row, var) pair of the LHS,
 *            3) d/dt = 0 for the algebraic variables (recomputed at each stage).
 *
 * @param u unknown vector
 * @param du_dt unkwon time derivative vector
 */
template <class T, int DIM>
void TransientOperator<T, DIM>::Mult(const mfem::Vector& u, mfem::Vector& du_dt) const {
  Catch_Time_Section("TransientOperator::Mult");

  const int fes_size = this->block_trueOffsets_.Size() - 1;

  // Working copy of the state (the algebraic variables of a split system may be projected)
  mfem::BlockVector ub(this->block_trueOffsets_);
  static_cast<mfem::Vector&>(ub) = u;

  // One vector per block (non-owning views, split along the block offsets)
  std::vector<mfem::Vector> v_vect;
  v_vect.reserve(static_cast<size_t>(fes_size));  // no reallocation -> views are not deep-copied
  for (int i = 0; i < fes_size; ++i) {
    mfem::Vector& bi = ub.GetBlock(i);
    v_vect.emplace_back(bi.GetData(), bi.Size());
  }

  // Todo(cci) : try to do different because of not satisfying
  auto* self = const_cast<TransientOperator<T, DIM>*>(this);
  self->SetExplicitTransientParameters(this->solve_time_step_, v_vect);

  mfem::BlockVector bb(this->block_trueOffsets_);
  mfem::BlockVector sol(this->block_trueOffsets_);
  sol = 0.;

  if (this->split_explicit_) {
    // 1) algebraic variables consistent with the differential ones at this stage
    //    (skipped if this exact state has already been projected, e.g. first stage of a step)
    if (!this->is_already_projected(ub)) {
      self->project_algebraic_variables(ub);
    }

    // 2) differential variables: lumped(row, var) d(var)/dt = -R_row
    // Source terms are part of the RHS nonlinear form (SourceTerm integrator)
    this->RHS->Mult(ub, bb);
    bb.Neg();
    for (size_t p = 0; p < this->evol_.size(); ++p) {
      const auto [r, c] = this->evol_[p];
      const mfem::Vector& in = bb.GetBlock(r);
      const mfem::Vector& m = this->lumped_[p];
      mfem::Vector& out = sol.GetBlock(c);
      for (int k = 0; k < out.Size(); ++k) {
        out[k] = in[k] / m[k];
      }
      // Essential dofs of the differential variable: no evolution
      const auto& ess = this->ess_tdof_list_[static_cast<size_t>(c)];
      for (int k = 0; k < ess.Size(); ++k) {
        out[ess[k]] = 0.;
      }
    }
    // 3) algebraic variables: d/dt = 0 (sol already set to zero)
  } else {
    // Source terms are part of the RHS nonlinear form (SourceTerm integrator)
    this->RHS->Mult(ub, bb);
    bb.Neg();

    for (int i = 0; i < fes_size; ++i) {
      apply_solver(this->M_solver_[static_cast<size_t>(i)], bb.GetBlock(i), sol.GetBlock(i));
    }
  }
  du_dt = sol;
}

/**
 * @brief  Solve the Backward-Euler equation: k = f(phi + dt*k, t), for the unknown k.
 *
 * @param dt current time step
 * @param u unknown vector
 * @param du_dt unkwon time derivative vector
 */
template <class T, int DIM>
void TransientOperator<T, DIM>::ImplicitSolve(const double dt, const mfem::Vector& u,
                                              mfem::Vector& du_dt) {
  Catch_Time_Section("TransientOperator::ImplicitSolve");

  const auto sc = this->height_;
  mfem::Vector v(u.GetData(), sc);
  mfem::Vector dv_dt(du_dt.GetData(), sc);
  const int fes_size = this->block_trueOffsets_.Size() - 1;

  // Apply BCs first: the forms below are then built (and their coefficients evaluated) with the
  // Dirichlet values of the current stage
  {
    MATools::MATrace::start();
    Catch_Time_Section("ImplicitSolve::ApplyBCs");

    // Get the unknown vectors of the auxiliary variables
    std::vector<mfem::Vector> auxvars_unk;
    for (const auto& auxvar_vec : this->auxvariables_) {
      for (const auto& auxvar : auxvar_vec->getVariables()) {
        auxvars_unk.emplace_back(auxvar.get_unknown());
      }
    }
    for (int i = 0; i < fes_size; ++i) {
      const int off = this->block_trueOffsets_[i];
      const int sz = this->block_trueOffsets_[i + 1] - off;
      mfem::Vector v_i(u.GetData() + off, sz);
      this->bcs_[static_cast<size_t>(i)]->SetDirichletBoundaryConditions(
          v_i, this->coefficients_[static_cast<size_t>(i)], auxvars_unk);
    }
    MATools::MATrace::stop("ApplyBCs");
  }

  MATools::MATrace::start();
  {
    Catch_Time_Section("ImplicitSolve::SetTransientParams");
    std::vector<mfem::Vector> v_vect;
    v_vect.reserve(static_cast<size_t>(fes_size));
    for (int i = 0; i < fes_size; ++i) {
      const int off = this->block_trueOffsets_[i];
      const int sz = this->block_trueOffsets_[i + 1] - off;
      v_vect.emplace_back(u.GetData() + off, sz);
    }
    this->SetTransientParameters(dt, v_vect);
    reduced_oper->SetParameters(dt, &v);
  }
  MATools::MATrace::stop("SetTransientParams");

  // Zero right-hand side: source terms are part of the RHS nonlinear form (SourceTerm
  // integrator)
  mfem::BlockVector zero_rhs(this->block_trueOffsets_);
  zero_rhs = 0.0;

  {
    MATools::MATrace::start();
    Catch_Time_Section("ImplicitSolve::CallMult");
#ifdef MFEM_USE_EXCEPTIONS
    // A failure inside the Newton solver (e.g. non-finite residual norm, MFEM_VERIFY) is a failure
    // of the time-step. Without MFEM_USE_EXCEPTIONS, MFEM aborts the computation.
    try {
      this->newton_solver_->Mult(zero_rhs, dv_dt);
    } catch (const mfem::ErrorException& error) {
      MATools::MATrace::stop("Solve");
      this->free_memory();
      throw TimeStepFailure(std::string("Nonlinear solver failed: ") + error.what());
    }
#else
    this->newton_solver_->Mult(zero_rhs, dv_dt);
#endif
    MATools::MATrace::stop("Solve");
  }

  // Same values on every rank: the Newton solver is built with MPI_COMM_WORLD (see
  // SolverNewton::create_solver), so its residual norms are reduced over all the ranks
  const bool converged = this->newton_solver_->GetConverged();
  const int iterations = this->newton_solver_->GetNumIterations();
  const double final_norm = this->newton_solver_->GetFinalNorm();
  this->newton_iterations_ = std::max(this->newton_iterations_, iterations);

  this->free_memory();

  if (!converged) {
    std::ostringstream msg;
    msg << "Nonlinear solver did not converge (" << iterations << " iterations, final norm "
        << final_norm << ", implicit stage step " << dt << ").";
    throw TimeStepFailure(msg.str());
  }
}

/**
 * @brief Overload the solver used to invert the mass matrix
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param SOLVER
 */
template <class T, int DIM>
void TransientOperator<T, DIM>::overload_mass_solver(VSolverType SOLVER) {
  this->mass_solver_ = SOLVER;
}

/**
 * @brief Overload the solver used to invert the mass matrix with its parameters
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param SOLVER
 * @param s_params
 */
template <class T, int DIM>
void TransientOperator<T, DIM>::overload_mass_solver(VSolverType SOLVER,
                                                     const Parameters& s_params) {
  this->mass_solver_ = SOLVER;
  this->mass_solver_params_ = s_params;
}

/**
 * @brief Overload the preconditioner for the solver used to invert the mass matrix
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param PRECOND
 */
template <class T, int DIM>
void TransientOperator<T, DIM>::overload_mass_preconditioner(VSolverType PRECOND) {
  this->mass_precond_ = PRECOND;
}

/**
 * @brief Overload the preconditioner for the solver used to invert the mass matrix and its
 * parameters
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param PRECOND
 * @param p_params
 */
template <class T, int DIM>
void TransientOperator<T, DIM>::overload_mass_preconditioner(VSolverType PRECOND,
                                                             const Parameters& p_params) {
  this->mass_precond_ = PRECOND;
  this->mass_precond_params_ = p_params;
}

/**
 * @brief Set the default solver and preconditioner associated with the mass matrix
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 *
 */
template <class T, int DIM>
void TransientOperator<T, DIM>::set_default_mass_solver() {
  auto s_params = Parameters(Parameter("description", "Default Mass Solver"));
  auto p_params = Parameters(Parameter("description", "Default Mass preconditioner"));

  this->mass_solver_ = HypreSolverType::HYPRE_GMRES;
  this->mass_solver_params_ = s_params;
  this->mass_precond_ = HyprePreconditionerType::HYPRE_ILU;
  this->mass_precond_params_ = p_params;
}

/**
 * @brief Build integrator for LHS
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 * @param integrator
 * @param vun
 * @param vauxn
 * @param all_params
 * @return SlothNLFormIntegrator<Variables<T, DIM>>*
 */
template <class T, int DIM>
SlothNLFormIntegrator<Variables<T, DIM>>* TransientOperator<T, DIM>::get_lhs_integrator(
    const double dt, const std::string integrator, const std::vector<mfem::ParGridFunction>& vun,
    const std::vector<mfem::ParGridFunction>& vauxn, const Parameters& all_params) {
  switch (Integrators::from(integrator)) {
    case Integrators::TimeDerivative: {
      return new TimeNLFormIntegrator<Variables<T, DIM>>(
          this->geometry_, dt, vun, vauxn, all_params, this->auxvariables_, this->coefficients_);
      break;
    }
    case Integrators::HeatTimeDerivative: {
      return new HeatTimeNLFormIntegrator<Variables<T, DIM>>(
          this->geometry_, dt, vun, vauxn, all_params, this->auxvariables_, this->coefficients_);
      break;
    }
    case Integrators::SplitTimeDerivative: {
      return new TimeCHNLFormIntegrator<Variables<T, DIM>>(
          this->geometry_, dt, vun, vauxn, all_params, this->auxvariables_, this->coefficients_);
      break;
    }
    default:
      mfem::mfem_error("LHS Integrators not found. Please check your data.");
      return nullptr;
  }
}

/**
 * @brief Free memory
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 *
 */
template <class T, int DIM>
void TransientOperator<T, DIM>::free_memory() {
  this->newton_solver_.reset();
  delete this->rhs_solver_;
  this->rhs_solver_ = nullptr;
  delete this->reduced_oper;
  this->reduced_oper = nullptr;
  delete this->RHS;
  this->RHS = nullptr;
  delete this->LHS;
  this->LHS = nullptr;

  this->free_mass_matrices();
}
/**
 * @brief Retrieve a coefficient by type and identifier.
 *
 * This function searches the list of coefficients associated with the
 * specified block (blk) and returns the coefficient that matches
 * the given glossary type, identifier and if given, the boundary id.
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 *
 * @param blk     Index of the block.
 * @param type    Type of coefficient.
 * @param id      Identifier of the coefficient.
 * @param bdr_id  Optional boundary id for boundary coefficients.
 *
 * @return An std::optional containing the matching Coefficient if found;
 *         std::nullopt otherwise.
 */
template <class T, int DIM>
std::optional<Coefficient> TransientOperator<T, DIM>::get_coefficient(const int blk,
                                                                      GlossaryType type,
                                                                      unsigned int id) {
  Coefficients coefficients = this->coefficients_[blk];

  for (unsigned int i = 0; i < coefficients.size(); i++) {
    auto coef = coefficients[i];
    if (coef.get_type() == type && coef.get_id() == id) {
      return coef;
    }
  }
  return std::nullopt;
}

/**
 * @brief Get the explicit time Coefficients object
 *        The mass matrix is constant (built once) only if both coefficients are scalar.
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 */
template <class T, int DIM>
void TransientOperator<T, DIM>::get_explicit_time_coefficients() {
  this->constant_mass_matrix_ = true;
  for (unsigned int k = 0; k < 2; ++k) {
    auto coef = this->get_coefficient(0, GlossaryType::ExplicitTime, k);

    if (!coef.has_value()) {
      coef = Coefficient(Glossary::Default, 1.0);
    } else {
      if (!(*coef).is_scalar() && !(*coef).is_explicit()) {
        mfem::mfem_error(
            ("Coefficient objects of type ExplicitTime for Explicit solver are either "
             "scalar or explicit. Default constant coefficient equal to one is used."));
      }
      if (!(*coef).is_scalar()) {
        this->constant_mass_matrix_ = false;
      }
    }

    this->explicit_time_coefficients_.add(*coef);
  }
}

/**
 * @brief Set the time of the operator and of its coefficients.
 *
 * Called by the MFEM ODE solvers before each evaluation: for a Runge-Kutta scheme, before each
 * stage, with the time of the stage (t + c_i dt). The coefficients depending on time (including
 * the Dirichlet values given by coefficients) are then evaluated at the time of the stage.
 *
 * @param t time
 */
template <class T, int DIM>
void TransientOperator<T, DIM>::SetTime(const mfem::real_t t) {
  mfem::TimeDependentOperator::SetTime(t);
  this->set_time_coefficients(t);
  this->explicit_time_coefficients_.set_time(t);
}

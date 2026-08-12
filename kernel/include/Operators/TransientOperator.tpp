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
#include <memory>
#include <string>
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
 * object
 *
 * @tparam T
 * @tparam DIM
 * @param spatials
 * @param ode
 * @param source_term_name
 */
template <class T, int DIM>
TransientOperator<T, DIM>::TransientOperator(

    std::vector<SpatialDiscretization<T, DIM>*> spatials,
    const std::vector<std::string>& rhs_integrators, TimeScheme::value ode,
    const std::string lhs_integrator, std::vector<AnalyticalFunctions<DIM>> source_term_name)
    : OperatorBase<T, DIM>(rhs_integrators, spatials, source_term_name),
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
 * @brief Construct a new TransientOperator object
 * @tparam T
 * @tparam DIM
 * @param spatial
 * @param params
 * @param vars
 * @param ode
 * @param auxvars
 * @param source_term_name
 */
template <class T, int DIM>
TransientOperator<T, DIM>::TransientOperator(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                                             const std::vector<std::string>& rhs_integrators,
                                             const Parameters& params, TimeScheme::value ode,
                                             const std::string lhs_integrator,
                                             std::vector<AnalyticalFunctions<DIM>> source_term_name)
    : OperatorBase<T, DIM>(rhs_integrators, spatials, params, source_term_name),
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
    default:
      mfem::mfem_error(
          "TimeDiscretization::set_ODE_solver: EulerImplicit, EulerExplicit, RungeKutta4 are "
          "available");
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

  OperatorBase<T, DIM>::initialize(initial_time, time_step, vars, auxvars);

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
  if (this->enable_amr_) {
    // Resynchronize height/width/z/block_trueOffsets_
    this->UpdateAfterMeshChange();
    // Resynchronize the OTHER height/width copy, inherited via mfem::TimeDependentOperator,
    // which UpdateAfterMeshChange() does not touch.
    this->mfem::TimeDependentOperator::height = this->OperatorBase<T, DIM>::height;
    this->mfem::TimeDependentOperator::width = this->OperatorBase<T, DIM>::width;
    // Resize k (the ODE solver's internal work vector)
    this->ode_solver_->Init(*this);
  }
  //// Constructing array of offsets
  const size_t unk_size = vect_unk.size();
  this->solve_time_step_ = current_time_step;
  //// Constructing BlockVector
  mfem::BlockVector block_unk(this->block_trueOffsets_);
  for (size_t i = 0; i < unk_size; i++) {
    auto& unk_i = *(vect_unk[i]);
    mfem::Vector& bb = block_unk.GetBlock(i);
    bb = unk_i;
  }
  //// Call ODE solver
  this->current_time_ = current_time;
  this->current_dt_ = current_time_step;
  this->ode_solver_->Step(block_unk, next_time, current_time_step);

  //// Updating vect_unk
  for (size_t i = 0; i < unk_size; i++) {
    auto& unk_i = *(vect_unk[i]);
    const mfem::Vector& bb = block_unk.GetBlock(i);
    unk_i = bb;
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief build the mass matrix
 *
 * @tparam T
 * @tparam DIM
 */
template <class T, int DIM>
void TransientOperator<T, DIM>::build_mass_matrix(const std::vector<mfem::Vector>& u_vect) {
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

  this->M_solver_.clear();
  for (auto* mat : this->Mmat_) {
    delete mat;
  }
  this->Mmat_.clear();
  for (unsigned int i = 0; i < u_vect.size(); i++) {
    if (M != nullptr) {
      delete M;
    }
    ////////////////
    // Mass matrix (constant)
    ////////////////
    M = new mfem::ParBilinearForm(this->fes_[i]);

    // Build Mfem Coefficient from Sloth Coefficient
    // Sloth coefficients are either scalar or explicit (variable + {auxiliary variables})
    mfem::ParGridFunction un(this->fes_[i]);
    un.SetFromTrueDofs(u_vect[i]);
    auto coef_a_exp = MfemCoefficient(0, this->explicit_time_coefficients_, un, vauxn);
    auto coef_b_exp = MfemCoefficient(1, this->explicit_time_coefficients_, un, vauxn);
    mfem::ProductCoefficient mass_coefficient = mfem::ProductCoefficient(coef_a_exp, coef_b_exp);

    M->AddDomainIntegrator(new mfem::LumpedIntegrator(new mfem::MassIntegrator(mass_coefficient)));

    M->Assemble(0);
    M->Finalize(0);

    this->Mmat_.emplace_back(M->ParallelAssemble());
    std::unique_ptr<mfem::HypreParMatrix> Me(
        this->Mmat_[i]->EliminateRowsCols(this->ess_tdof_list_[i]));

    this->mass_matrix_solver_ =
        std::make_shared<LSolver>(this->mass_solver_, this->mass_solver_params_,
                                  this->mass_precond_, this->mass_precond_params_, *this->Mmat_[i]);
    this->M_solver_.emplace_back(this->mass_matrix_solver_->get_solver());
  }
}

/**
 * @brief Build the nonlinear form associated with LHS of the PDEs
 *
 * @tparam T
 * @tparam DIM
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
 * @tparam T
 * @tparam DIM
 * @param dt
 * @param u_vect
 */
template <class T, int DIM>
SlothNLFormIntegrator<Variables<T, DIM>>* TransientOperator<T, DIM>::set_lhs_nlfi_ptr(
    const double dt, const std::vector<mfem::Vector>& u) {
  Catch_Time_Section("TransientOperator::set_lhs_nlfi_ptr");

  std::vector<mfem::ParGridFunction> vun;
  std::vector<mfem::ParGridFunction> vauxn;
  for (unsigned int i = 0; i < u.size(); i++) {
    mfem::ParGridFunction un(this->fes_[i]);
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
  lhs->init();
  return lhs;
}

/**
 * @brief Set current dt, unk values - needed to compute action and Jacobian.
 *solution_coef
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
 * @tparam T
 * @tparam DIM
 */
template <class T, int DIM>
void TransientOperator<T, DIM>::SetExplicitTransientParameters(
    const double dt, const std::vector<mfem::Vector>& un_vect) {
  Catch_Time_Section("TransientOperator::SetExplicitTransientParameters");
  ////////////////////////////////////////////
  // Variable mass matrix
  ////////////////////////////////////////////
  this->build_mass_matrix(un_vect);

  ////////////////////////////////////////////
  // PhaseField non linear form
  ////////////////////////////////////////////
  this->build_rhs_nonlinear_form(dt, un_vect);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/**
 * @brief Compute the right-hand side of the ODE system.
 *
 * @param u unknown vector
 * @param du_dt unkwon time derivative vector
 */
template <class T, int DIM>
void TransientOperator<T, DIM>::Mult(const mfem::Vector& u, mfem::Vector& du_dt) const {
  Catch_Time_Section("TransientOperator::Mult");

  const auto sc = this->height_;
  mfem::Vector v(u.GetData(), sc);
  mfem::Vector dv_dt(du_dt.GetData(), sc);

  // Todo(cci) change with BlockVector
  std::vector<mfem::Vector> v_vect;
  v_vect.emplace_back(v);

  // Todo(cci) : try to do different because of not satisfying
  const_cast<TransientOperator<T, DIM>*>(this)->SetExplicitTransientParameters(
      this->solve_time_step_, v_vect);

  const int fes_size = this->block_trueOffsets_.Size() - 1;
  mfem::BlockVector bb(this->block_trueOffsets_);

  this->RHS->Mult(v, bb);

  // Source term
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
    bb -= source_term;
  }
  bb.Neg();

  mfem::BlockVector sol(this->block_trueOffsets_);
  sol = 0.;

  for (int i = 0; i < fes_size; ++i) {
    mfem::Vector& soli = sol.GetBlock(i);
    mfem::Vector& bb_i = bb.GetBlock(i);
    std::visit(
        [&](auto&& arg) {
          using TT = std::decay_t<decltype(arg)>;
          if constexpr (!std::is_same_v<TT, std::shared_ptr<std::monostate>>) {
            if (bb_i.Size() != arg->NumCols()) {
              std::string msg = "Size mismatch: bb_i.Size() [" + std::to_string(bb_i.Size()) +
                                "] != arg->NumCols() [" + std::to_string(arg->NumCols()) + "]";
              MFEM_ABORT(msg.c_str());
            }
            arg->Mult(bb_i, soli);
          }
        },
        this->M_solver_[i]);
  }
  dv_dt = sol;
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

  MATools::MATrace::start();
  {
    Catch_Time_Section("ImplicitSolve::SetTransientParams");
    std::vector<mfem::Vector> v_vect;
    auto sc_1 = 0;
    auto sc_2 = sc / fes_size;
    for (int i = 0; i < fes_size; ++i) {
      mfem::Vector v_i(u.GetData() + sc_1, sc_2);
      sc_1 += sc_2;
      v_vect.emplace_back(v_i);
    }
    this->SetTransientParameters(dt, v_vect);
  }
  MATools::MATrace::stop("SetTransientParams");

  // Apply BCs
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
    auto sc_1 = 0;
    auto sc_2 = sc / fes_size;
    for (int i = 0; i < fes_size; ++i) {
      mfem::Vector v_i(u.GetData() + sc_1, sc_2);
      this->bcs_[i]->SetDirichletBoundaryConditions(v_i, this->coefficients_[i], auxvars_unk);
      sc_1 += sc_2;
    }
    reduced_oper->SetParameters(dt, &v);
    MATools::MATrace::stop("ApplyBCs");
  }

  // Source term
  mfem::BlockVector source_term(this->block_trueOffsets_);
  source_term = 0.0;
  {
    MATools::MATrace::start();
    Catch_Time_Section("ImplicitSolve::SourceTerm");
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
    MATools::MATrace::stop("SourceTerm");
  }

  {
    MATools::MATrace::start();
    Catch_Time_Section("ImplicitSolve::CallMult");
    this->newton_solver_->Mult(source_term, dv_dt);
    MATools::MATrace::stop("Solve");
  }

  bool converged = this->newton_solver_->GetConverged();

  this->free_memory();

  MFEM_VERIFY(converged, "Nonlinear solver did not converge.");
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

  for (auto* mat : this->Mmat_) {
    delete mat;
  }
  this->Mmat_.clear();
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
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 */
template <class T, int DIM>
void TransientOperator<T, DIM>::get_explicit_time_coefficients() {
  for (int k = 0; k < 2; ++k) {
    auto coef = this->get_coefficient(0, GlossaryType::ExplicitTime, k);

    if (!coef.has_value()) {
      coef = Coefficient(Glossary::Default, 1.0);
    } else {
      if (!(*coef).is_scalar() && !(*coef).is_explicit()) {
        mfem::mfem_error(
            ("Coefficient objects of type ExplicitTime for Explicit solver are either "
             "scalar or explicit. Default constant coefficient equal to one is used."));
      }
    }

    this->explicit_time_coefficients_.add(*coef);
  }
}

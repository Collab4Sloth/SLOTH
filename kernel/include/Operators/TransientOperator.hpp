/**
 * @file TransientOperator.hpp
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

#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "AnalyticalFunctions/AnalyticalFunctions.hpp"
#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Operators/OperatorBase.hpp"
#include "Operators/ReducedOperator.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Solvers/SlothBDF2Solver.hpp"
#include "Solvers/LSolver.hpp"
#include "Spatial/Spatial.hpp"
#include "Utils/Utils.hpp"
#include "Variables/Variables.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief TransientOperator class
 *
 */
template <class T, int DIM>
class TransientOperator : public OperatorBase<T, DIM>, public mfem::TimeDependentOperator {
 private:
  double solve_time_step_{std::numeric_limits<double>::quiet_NaN()};
  std::shared_ptr<mfem::ODESolver> ode_solver_;
  // BDF2 scheme (nullptr otherwise): history and lagged auxiliary variables handled by the solver
  std::shared_ptr<SlothBDF2Solver<T, DIM>> bdf2_solver_;
  // BDF2: previous time step saved at the beginning of a time-step (save_time_scheme_state)
  mfem::real_t saved_previous_time_step_{0.};
  // Other schemes: lagged auxiliary variables of a partitioned coupling reported once
  Variables<T, DIM>* vars_{nullptr};
  bool lagged_aux_info_printed_{false};
  void report_lagged_auxiliary_variables();

  bool is_explicit_;
  void set_ODE_solver(const TimeScheme::value& ode_solver);

  VSolverType mass_solver_;
  VSolverType mass_precond_;
  Parameters mass_solver_params_;
  Parameters mass_precond_params_;

  void set_default_mass_solver();
  SlothNLFormIntegrator<Variables<T, DIM>>* set_lhs_nlfi_ptr(const double dt,
                                                             const std::vector<mfem::Vector>& u);

  SlothNLFormIntegrator<Variables<T, DIM>>* get_lhs_integrator(
      const double dt, const std::string integrator, const std::vector<mfem::ParGridFunction>& vun,
      const std::vector<mfem::ParGridFunction>& vauxn, const Parameters& all_params);
  std::string lhs_integrator_;

  void free_memory();
  void free_mass_matrices();

  Coefficients explicit_time_coefficients_;
  std::optional<Coefficient> get_coefficient(const int blk, GlossaryType type, unsigned int id);
  void get_explicit_time_coefficients();

  //////////////////////////////////////////////////////////
  // Explicit time scheme for DAE system
  //////////////////////////////////////////////////////////
  // Split systems (index-1 DAE, e.g. split Cahn-Hilliard / Allen-Cahn): some rows carry no
  // time derivative (algebraic constraints). With an explicit scheme ("Split" in the name of
  // the LHS integrator), the system is solved in state-space form:
  //   - the LHS integrator declares which row carries d(variable)/dt (time_derivative_map());
  //   - at each stage, the algebraic variables are projected onto the constraints (Newton,
  //     Jacobian kept between projections), then only the differential variables are advanced
  //     (lumped mass of the LHS), with d/dt = 0 for the algebraic variables;
  //   - the algebraic variables are projected again at the end of each step (consistent output).
  // Requirements: constraint Jacobian invertible w.r.t. the algebraic variables (index 1), one
  // constraint per algebraic variable, no Dirichlet condition on algebraic variables, RHS
  // integrators reading the algebraic variables from the current state (elfun).
  //////////////////////////////////////////////////////////

  // Apply a SLOTH linear solver (variant) to a vector
  static void apply_solver(const VSharedMFEMSolver& solver, const mfem::Vector& in,
                           mfem::Vector& out);

  /// Split (DAE) systems with explicit schemes ("Split" in the name of the LHS integrator)
  // true once the structure of the system has been analysed
  bool split_checked_{false};
  // true if explicit scheme and split LHS integrator
  bool split_explicit_{false};
  // (row, variable): row 'row' carries d(variable)/dt (given by the LHS integrator)
  std::vector<std::pair<int, int>> evol_;
  // The constraint Jacobian is rebuilt if a correction reduces the residual by less than this
  // factor (optional parameter "constraint_jacobian_min_reduction", default 0.1)
  double constraint_min_reduction_{0.1};
  // (row, variable): constraint row 'row' is solved for the algebraic variable 'variable'
  std::vector<std::pair<int, int>> alg_;
  // lumped mass (row sums of the LHS Jacobian block (row, variable)) of each evol_ pair
  std::vector<mfem::Vector> lumped_;
  // true if every ExplicitTime coefficient is scalar (or absent): lumped masses built once
  bool split_mass_constant_{false};
  // Signature (state, auxiliary variables, time, time step) of the last projected state
  mfem::Vector last_projection_signature_;
  // Constraint Jacobian blocks (copies) and their solvers, kept between projections
  // (modified Newton: refreshed only if the convergence is too slow, or after a mesh change).
  // Declared matrices first: solvers are destroyed before the matrices they refer to.
  std::vector<std::shared_ptr<mfem::HypreParMatrix>> constraint_mats_;
  std::vector<std::shared_ptr<LSolver>> constraint_solvers_;
  void build_constraint_solvers(const mfem::BlockVector& ub);
  void free_constraint_solvers();
  void setup_split_structure(const double dt, const std::vector<mfem::Vector>& u_vect);
  void compute_split_lumped_masses();
  void project_algebraic_variables(mfem::BlockVector& ub);
  mfem::Vector projection_signature(const mfem::Vector& u) const;
  bool is_already_projected(const mfem::Vector& u) const;

  //////////////////////////////////////////////////////////
 protected:
  /// Mass operator

  mfem::ParBilinearForm* M;  // mass operator
  std::shared_ptr<LSolver> mass_matrix_solver_;
  std::vector<std::shared_ptr<LSolver>> mass_solvers_;  // one mass solver per block
  std::vector<VSharedMFEMSolver> M_solver_;  // Krylov solver for inverting the mass matrix M
  /// Left-Hand-Side
  mfem::ParBlockNonlinearForm* LHS;

  std::vector<mfem::HypreParMatrix*> Mmat_;

  void build_mass_matrix(const std::vector<mfem::Vector>& u_vect);
  void build_lhs_nonlinear_form(const double dt, const std::vector<mfem::Vector>& u);
  // true if both ExplicitTime coefficients are scalar: mass matrix built once (explicit schemes)
  bool constant_mass_matrix_{true};

  // Reduced Operator
  PhaseFieldReducedOperator* reduced_oper;

 public:
  // Time scheme (adaptive time-stepping, see OperatorBase)
  double get_max_time_step_ratio() const;
  void save_time_scheme_state();
  void restore_time_scheme_state();

  TransientOperator(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                    const std::vector<std::string>& rhs_integrators, TimeScheme::value ode,
                    const std::string lhs_integrator);

  TransientOperator(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                    const std::vector<std::string>& rhs_integrators, const Parameters& params,
                    TimeScheme::value ode, const std::string lhs_integrator);

  ~TransientOperator() override;

  void Mult(const mfem::Vector& u, mfem::Vector& du_dt) const override;
  // Called by the MFEM ODE solvers before each evaluation (each stage of a Runge-Kutta
  // scheme): the coefficients are set to the time of the stage.
  void SetTime(const mfem::real_t t) override;
  void ImplicitSolve(const double dt, const mfem::Vector& u, mfem::Vector& k) override;

  // User-defined Solvers
  void overload_mass_solver(VSolverType SOLVER);
  void overload_mass_solver(VSolverType SOLVER, const Parameters& s_params);
  void overload_mass_preconditioner(VSolverType PRECOND);
  void overload_mass_preconditioner(VSolverType PRECOND, const Parameters& p_params);

  void SetExplicitTransientParameters(const double dt, const std::vector<mfem::Vector>& un_vect);

  // Virtual methods
  void initialize(const double& initial_time, const double time_step, Variables<T, DIM>& vars,
                  std::vector<Variables<T, DIM>*> auxvars) override;
  void SetTransientParameters(const double dt, const std::vector<mfem::Vector>& u_vect) override;
  void solve(std::vector<std::unique_ptr<mfem::Vector>>& vect_unk, double& next_time,
             const double& current_time, double current_time_step, const int iter) override;
};

#include "Operators/TransientOperator.tpp"
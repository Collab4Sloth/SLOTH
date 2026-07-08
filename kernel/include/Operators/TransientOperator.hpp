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

#include <memory>
#include <string>
#include <vector>

#include "AnalyticalFunctions/AnalyticalFunctions.hpp"
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

/**
 * @brief TransientOperator class
 *
 */
template <class T, int DIM>
class TransientOperator : public OperatorBase<T, DIM>, public mfem::TimeDependentOperator {
 private:
  double solve_time_step_{std::numeric_limits<double>::quiet_NaN()};
  std::shared_ptr<mfem::ODESolver> ode_solver_;

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

  Coefficients explicit_time_coefficients_;
  std::optional<Coefficient> get_coefficient(const int blk, GlossaryType type, unsigned int id);
  void get_explicit_time_coefficients();

 protected:
  /// Mass operator

  mfem::ParBilinearForm* M;  // mass operator
  std::shared_ptr<LSolver> mass_matrix_solver_;
  std::vector<VSharedMFEMSolver> M_solver_;  // Krylov solver for inverting the mass matrix )M
  /// Left-Hand-Side
  mfem::ParBlockNonlinearForm* LHS;

  std::vector<mfem::HypreParMatrix*> Mmat_;

  void build_mass_matrix(const std::vector<mfem::Vector>& u_vect);
  void build_lhs_nonlinear_form(const double dt, const std::vector<mfem::Vector>& u);
  bool constant_mass_matrix_{true};

  // Reduced Operator
  PhaseFieldReducedOperator* reduced_oper;

 public:
  TransientOperator(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                    const std::vector<std::string>& rhs_integrators, TimeScheme::value ode,
                    const std::string lhs_integrator);

  TransientOperator(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                    const std::vector<std::string>& rhs_integrators, TimeScheme::value ode,
                    const std::string lhs_integrator,
                    std::vector<AnalyticalFunctions<DIM>> source_term_name);

  TransientOperator(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                    const std::vector<std::string>& rhs_integrators, const Parameters& params,
                    TimeScheme::value ode, const std::string lhs_integrator);

  TransientOperator(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                    const std::vector<std::string>& rhs_integrators, const Parameters& params,
                    TimeScheme::value ode, const std::string lhs_integrator,
                    std::vector<AnalyticalFunctions<DIM>> source_term_name);

  void Mult(const mfem::Vector& u, mfem::Vector& du_dt) const override;
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

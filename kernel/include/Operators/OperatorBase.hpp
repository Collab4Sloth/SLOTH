/**
 * @file OperatorBase.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Base class for building Steady and TimeDependent PhaseFieldOperators
 * @version 0.1
 * @date 2025-09-05
 *
 * @copyright CEA (C) 2025
 *
 * @anchor OperatorBase
 *
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
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "AnalyticalFunctions/AnalyticalFunctions.hpp"
#include "BCs/BoundaryConditions.hpp"
#include "Coefficients/Coefficients.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Solvers/LSolver.hpp"
#include "Solvers/NLSolver.hpp"
#include "Spatial/Spatial.hpp"
#include "Utils/Utils.hpp"
#include "Variables/Variables.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Base class for building Steady and TimeDependent PhaseFieldOperators
 *
 * @tparam T Finite Element collection (mfem object)
 * @tparam DIM Spatial dimension
 */
template <class T, int DIM>
class OperatorBase : public mfem::Operator {
 private:
  T* fecollection_;

  NLSolverType nl_solver_;
  VSolverType solver_;
  VSolverType precond_;
  Parameters solver_params_;
  Parameters precond_params_;
  void set_default_solver();

 protected:
  Geometry geometry_ = Geometry::Cartesian;
  Parameters nl_solver_params_;
  std::optional<Coefficient> energy_coefficient_;
  std::optional<Coefficient> grad_energy_coefficient_;
  std::vector<Coefficients> coefficients_;
  std::vector<Variables<T, DIM>*> auxvariables_;
  std::string description_{"UNKNOWN OPERATOR"};
  const Parameter default_p_ = Parameter("default parameter", false);
  const Parameters default_params_ = Parameters(default_p_);
  const Parameters& params_;
  /// Time integral results
  std::multimap<IterationKey, SpecializedValue> time_specialized_;
  std::map<std::string, std::multimap<IterationKey, SpecializedValue>> time_iso_specialized_;
  /// Finite Element Space and BCs
  mfem::Array<mfem::ParFiniteElementSpace*> fes_;
  std::vector<mfem::Array<int>> ess_tdof_list_;
  mfem::Array<int> block_trueOffsets_;

  /// Right-Hand-Side
  mfem::ParBlockNonlinearForm* RHS;

  NLSolver* rhs_solver_;
  std::shared_ptr<mfem::NewtonSolver> newton_solver_;

  /// Boundary conditions
  std::vector<BoundaryConditions<T, DIM>*> bcs_;
  std::vector<mfem::Array<int>> array_bdr_;

  std::vector<std::function<double(const mfem::Vector&, double)>> src_func_;

  double current_dt_;
  double current_time_;
  int height_;
  mutable mfem::Vector z;  // auxiliary vector

  std::vector<std::string> rhs_integrators_;

  void build_rhs_nonlinear_form(const std::vector<mfem::Vector>& u);
  void SetNewtonAlgorithm(mfem::Operator* oper);

  int compute_total_height(const std::vector<SpatialDiscretization<T, DIM>*>& spatials);
  int compute_total_width(const std::vector<SpatialDiscretization<T, DIM>*>& spatials);
  SlothNLFormIntegrator<Variables<T, DIM>>* get_rhs_integrator(
      const std::string integrator, const std::vector<mfem::ParGridFunction>& vun,
      const std::vector<mfem::ParGridFunction>& vauxn, const Parameters& all_params);
  SlothNLFormIntegrator<Variables<T, DIM>>* get_bdr_integrator(
      const std::string integrator, const std::vector<mfem::ParGridFunction>& vun,
      const std::vector<mfem::ParGridFunction>& vauxn, const Parameters& all_params,
      const unsigned int block, const unsigned int bdr_id);

 public:
  OperatorBase(const std::vector<std::string>& integrators,
               std::vector<SpatialDiscretization<T, DIM>*> spatials);

  OperatorBase(const std::vector<std::string>& integrators,
               std::vector<SpatialDiscretization<T, DIM>*> spatials,
               const std::vector<AnalyticalFunctions<DIM>>& source_term_name);
  OperatorBase(const std::vector<std::string>& integrators,
               std::vector<SpatialDiscretization<T, DIM>*> spatials, const Parameters& params);

  OperatorBase(const std::vector<std::string>& integrators,
               std::vector<SpatialDiscretization<T, DIM>*> spatials, const Parameters& params,
               const std::vector<AnalyticalFunctions<DIM>>& source_term_name);

  void set_coefficients(const std::vector<Coefficients>& coefficients,
                        bool enable_compute_energies);

  inline double get_current_time_step() const { return current_dt_; }

  void ComputeError(const int& it, const double& t, const double& dt, const int id_var,
                    const std::string& name, const mfem::Vector& u,
                    std::function<double(const mfem::Vector&, double)> solution_func);
  void ComputeIntegral(const int& it, const double& t, const double& dt, const int id_var,
                       const std::string& name, const mfem::Vector& u, const double lower_bound,
                       const double upper_bound);
  void ComputeEnergies(const int& it, const double& t, const double& dt,
                       const std::vector<mfem::Vector>& u);

  void ComputeIsoVal(const int& it, const double& t, const double& dt, const int var_id,
                     const std::string& var_name, const mfem::Vector& u, const double& iso_value);
  void get_source_term(const int id_block,
                       const std::function<double(const mfem::Vector&, double)>& src_func,
                       mfem::Vector& source_term, mfem::ParLinearForm* RHHS) const;

  const std::multimap<IterationKey, SpecializedValue>& get_time_specialized() const;
  const std::map<std::string, std::multimap<IterationKey, SpecializedValue>>
  get_time_iso_specialized() const;

  void clear_time_specialized();
  void clear_iso_time_specialized();

  virtual ~OperatorBase() = default;

  std::string get_description() { return this->description_; }

  // User-defined Solvers
  void overload_nl_solver(NLSolverType NLSOLVER);
  void overload_nl_solver(NLSolverType NLSOLVER, const Parameters& nl_params);

  void overload_solver(VSolverType SOLVER);
  void overload_solver(VSolverType SOLVER, const Parameters& s_params);

  void overload_preconditioner(VSolverType PRECOND);
  void overload_preconditioner(VSolverType PRECOND, const Parameters& p_params);

  void setGeometry(Geometry geometry);
  // Virtual methods
  virtual void initialize(const double& initial_time, Variables<T, DIM>& vars,
                          std::vector<Variables<T, DIM>*> auxvars);

  // Pure virtual methods
  virtual void SetTransientParameters(const std::vector<mfem::Vector>& u_vect) = 0;
  virtual void solve(std::vector<std::unique_ptr<mfem::Vector>>& vect_unk, double& next_time,
                     const double& current_time, double current_time_step, const int iter) = 0;
  SlothNLFormIntegrator<Variables<T, DIM>>* set_nlfi_ptr(const std::string nlfi,
                                                         const std::vector<mfem::Vector>& u);
  SlothNLFormIntegrator<Variables<T, DIM>>* set_bdr_nlfi_ptr(const std::string nlfi,
                                                             const std::vector<mfem::Vector>& u,
                                                             const unsigned int block,
                                                             const unsigned int bdr_id);
  void set_time_coefficients(double time);
};
#include "Operators/OperatorBase.tpp"

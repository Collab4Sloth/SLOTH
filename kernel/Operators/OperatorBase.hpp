/**
 * @file OperatorBase.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Base class for building Steady and TimeDependent PhaseFieldOperators
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
#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "AnalyticalFunctions/AnalyticalFunctions.hpp"
#include "BCs/BoundaryConditions.hpp"
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

#pragma once

/**
 * @brief Base class for building Steady and TimeDependent PhaseFieldOperators
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
class OperatorBase : public mfem::Operator {
 private:
  T* fecollection_;

  NLSolverType nl_solver_;
  VSolverType solver_;
  VSolverType precond_;
  Parameters nl_solver_params_;
  Parameters solver_params_;
  Parameters precond_params_;
  void set_default_solver();

 protected:
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

  std::vector<std::function<double(const mfem::Vector&, double)>> src_func_;

  double current_dt_;
  double current_time_;
  int height_;
  mutable mfem::Vector z;  // auxiliary vector
  NLFI* nlfi_ptr_;

  void build_rhs_nonlinear_form(const double dt, const std::vector<mfem::Vector>& u);
  void SetNewtonAlgorithm(mfem::Operator* oper);

  int compute_total_height(const std::vector<SpatialDiscretization<T, DIM>*>& spatials);
  int compute_total_width(const std::vector<SpatialDiscretization<T, DIM>*>& spatials);

 public:
  explicit OperatorBase(std::vector<SpatialDiscretization<T, DIM>*> spatials);

  OperatorBase(std::vector<SpatialDiscretization<T, DIM>*> spatials,
               const std::vector<AnalyticalFunctions<DIM>>& source_term_name);
  OperatorBase(std::vector<SpatialDiscretization<T, DIM>*> spatials, const Parameters& params);

  OperatorBase(std::vector<SpatialDiscretization<T, DIM>*> spatials, const Parameters& params,
               const std::vector<AnalyticalFunctions<DIM>>& source_term_name);

  void ComputeError(const int& it, const double& t, const double& dt, const int id_var,
                    const std::string& name, const mfem::Vector& u,
                    std::function<double(const mfem::Vector&, double)> solution_func);
  void ComputeIntegral(const int& it, const double& t, const double& dt, const int id_var,
                       const std::string& name, const mfem::Vector& u, const double lower_bound,
                       const double upper_bound);
  void ComputeIsoVal(const int& it, const double& t, const double& dt, const int var_id,
                     const std::string& var_name, const mfem::Vector& u, const double& iso_value);
  void get_source_term(const int id_block,
                       const std::function<double(const mfem::Vector&, double)>& src_func,
                       mfem::Vector& source_term, mfem::ParLinearForm* RHHS) const;

  const std::multimap<IterationKey, SpecializedValue> get_time_specialized() const;
  const std::map<std::string, std::multimap<IterationKey, SpecializedValue>>
  get_time_iso_specialized() const;

  void clear_time_specialized();
  void clear_iso_time_specialized();

  virtual ~OperatorBase();

  std::string get_description() { return this->description_; }

  // User-defined Solvers
  void overload_nl_solver(NLSolverType NLSOLVER);
  void overload_nl_solver(NLSolverType NLSOLVER, const Parameters& nl_params);

  void overload_solver(VSolverType SOLVER);
  void overload_solver(VSolverType SOLVER, const Parameters& s_params);

  void overload_preconditioner(VSolverType PRECOND);
  void overload_preconditioner(VSolverType PRECOND, const Parameters& p_params);

  // Virtual methods
  // virtual void initialize(const double &initial_time, Variables<T, DIM> &vars);
  virtual void initialize(const double& initial_time, Variables<T, DIM>& vars,
                          std::vector<Variables<T, DIM>*> auxvars);

  // Pure virtual methods
  virtual void set_default_properties() = 0;
  virtual void SetConstantParameters(const double dt, const std::vector<mfem::Vector>& u_vect) = 0;
  virtual void SetTransientParameters(const double dt, const std::vector<mfem::Vector>& u_vect) = 0;
  virtual void solve(std::vector<std::unique_ptr<mfem::Vector>>& vect_unk, double& next_time,
                     const double& current_time, double current_time_step, const int iter) = 0;
  virtual NLFI* set_nlfi_ptr(const double dt, const std::vector<mfem::Vector>& u) = 0;
  virtual void get_parameters() = 0;
  virtual void ComputeEnergies(const int& it, const double& dt, const double& t,
                               const std::vector<mfem::Vector>& u) = 0;
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Return the total height (output=rows of Operator) of the PDE system
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param spatials
 * @return int
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
int OperatorBase<T, DIM, NLFI, LHS_NLFI>::compute_total_height(
    const std::vector<SpatialDiscretization<T, DIM>*>& spatials) {
  int total_size = 0;
  for (const auto* s : spatials) {
    total_size += s->getSize();
  }
  return total_size;
}

/**
 * @brief Return the total width (input=column of Operator) of the PDE system
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @tparam LHS_NLFI
 * @param spatials
 * @return int
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
int OperatorBase<T, DIM, NLFI, LHS_NLFI>::compute_total_width(
    const std::vector<SpatialDiscretization<T, DIM>*>& spatials) {
  // return spatials.size();
  int total_size = 0;
  for (const auto* s : spatials) {
    total_size += s->getSize();
  }
  return total_size;
}

/**
 * @brief Construct a new OperatorBase object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param spatial
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
OperatorBase<T, DIM, NLFI, LHS_NLFI>::OperatorBase(
    std::vector<SpatialDiscretization<T, DIM>*> spatials)
    : mfem::Operator(this->compute_total_height(spatials), this->compute_total_width(spatials)),
      params_(default_params_),
      RHS(NULL),
      current_dt_(0.0),
      current_time_(0.0),
      height_(height),
      z(height),
      nlfi_ptr_(nullptr) {
  this->fes_.SetSize(spatials.size());
  this->bcs_.reserve(spatials.size());
  this->ess_tdof_list_.reserve(spatials.size());
  this->block_trueOffsets_.SetSize(spatials.size() + 1);

  this->block_trueOffsets_[0] = 0;
  int i = 0;
  for (const auto* s : spatials) {
    this->fes_[i] = s->get_finite_element_space();
    this->block_trueOffsets_[i + 1] = this->fes_[i]->GetTrueVSize();
    i++;
  }
  this->block_trueOffsets_.PartialSum();
  this->set_default_solver();
}

/**
 * @brief Construct a new OperatorBase object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param spatial
 * @param source_term_name
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
OperatorBase<T, DIM, NLFI, LHS_NLFI>::OperatorBase(
    std::vector<SpatialDiscretization<T, DIM>*> spatials,
    const std::vector<AnalyticalFunctions<DIM>>& source_term_name)
    : mfem::Operator(this->compute_total_height(spatials), this->compute_total_width(spatials)),
      params_(default_params_),
      RHS(NULL),
      current_dt_(0.0),
      current_time_(0.0),
      height_(height),
      z(height),
      nlfi_ptr_(nullptr) {
  this->fes_.SetSize(spatials.size());
  this->bcs_.reserve(spatials.size());
  this->ess_tdof_list_.reserve(spatials.size());

  this->block_trueOffsets_.SetSize(spatials.size() + 1);

  this->block_trueOffsets_[0] = 0;
  int i = 0;
  for (const auto* s : spatials) {
    this->fes_[i] = s->get_finite_element_space();
    this->block_trueOffsets_[i + 1] = this->fes_[i]->GetTrueVSize();
    i++;
  }
  this->block_trueOffsets_.PartialSum();
  for (const auto& src : source_term_name) {
    this->src_func_.emplace_back(std::move(src.getFunction()));
  }
  this->set_default_solver();
}

/**
 * @brief Construct a new OperatorBase object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param spatial
 * @param params
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
OperatorBase<T, DIM, NLFI, LHS_NLFI>::OperatorBase(
    std::vector<SpatialDiscretization<T, DIM>*> spatials, const Parameters& params)
    : mfem::Operator(this->compute_total_height(spatials), this->compute_total_width(spatials)),
      params_(params),
      RHS(NULL),
      current_dt_(0.0),
      current_time_(0.0),
      height_(height),
      z(height),
      nlfi_ptr_(nullptr) {
  this->fes_.SetSize(spatials.size());
  this->bcs_.reserve(spatials.size());
  this->ess_tdof_list_.reserve(spatials.size());

  this->block_trueOffsets_.SetSize(spatials.size() + 1);

  this->block_trueOffsets_[0] = 0;
  int i = 0;
  for (const auto* s : spatials) {
    this->fes_[i] = s->get_finite_element_space();
    this->block_trueOffsets_[i + 1] = this->fes_[i]->GetTrueVSize();
    i++;
  }
  this->block_trueOffsets_.PartialSum();
  this->set_default_solver();
}

/**
 * @brief  Construct a new Phase Field Operator Base< T,  DIM, NLFI>:: Phase Field Operator
 * Base object
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param spatial
 * @param params
 * @param vars
 * @param auxvars
 * @param source_term_name
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
OperatorBase<T, DIM, NLFI, LHS_NLFI>::OperatorBase(
    std::vector<SpatialDiscretization<T, DIM>*> spatials, const Parameters& params,
    const std::vector<AnalyticalFunctions<DIM>>& source_term_name)
    : mfem::Operator(this->compute_total_height(spatials), this->compute_total_width(spatials)),
      params_(params),
      RHS(NULL),
      current_dt_(0.0),
      current_time_(0.0),
      height_(height),
      z(height),
      nlfi_ptr_(nullptr) {
  this->fes_.SetSize(spatials.size());
  this->bcs_.reserve(spatials.size());
  this->ess_tdof_list_.reserve(spatials.size());

  this->block_trueOffsets_.SetSize(spatials.size() + 1);

  this->block_trueOffsets_[0] = 0;
  int i = 0;
  for (const auto* s : spatials) {
    this->fes_[i] = s->get_finite_element_space();
    this->block_trueOffsets_[i + 1] = this->fes_[i]->GetTrueVSize();
    i++;
  }
  this->block_trueOffsets_.PartialSum();
  for (const auto& src : source_term_name) {
    auto s = src.getFunction();
    this->src_func_.emplace_back(s);
  }
  this->set_default_solver();
}

/**
 * @brief  Initialization stage (call by imeDiscretization<PST, OPE, VAR>::initialize())
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param initial_time
 * @param vars
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void OperatorBase<T, DIM, NLFI, LHS_NLFI>::initialize(const double& initial_time,
                                                      Variables<T, DIM>& vars,
                                                      std::vector<Variables<T, DIM>*> auxvars) {
  Catch_Time_Section("OperatorBase::initialize");

  this->auxvariables_ = auxvars;
  const auto nvars = vars.get_variables_number();
  std::vector<mfem::Vector> u_vect;
  u_vect.reserve(nvars);
  for (auto iv = 0; iv < nvars; iv++) {
    auto& vv = vars.getIVariable(iv);
    auto u = vv.get_unknown();

    this->bcs_.emplace_back(vv.get_boundary_conditions());
    this->ess_tdof_list_.emplace_back(this->bcs_[iv]->GetEssentialDofs());
    this->bcs_[iv]->SetBoundaryConditions(u);
    vv.update(u);
    u_vect.emplace_back(u);
  }
  this->SetConstantParameters(this->current_dt_, u_vect);
  this->SetTransientParameters(this->current_dt_, u_vect);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief Build the NonLinear Form Integrator associated with the RHS of the PDEs
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param dt
 * @param u
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void OperatorBase<T, DIM, NLFI, LHS_NLFI>::build_rhs_nonlinear_form(
    const double dt, const std::vector<mfem::Vector>& u_vect) {
  if (this->RHS != nullptr) {
    delete this->RHS;
  }
  this->RHS = new mfem::ParBlockNonlinearForm(this->fes_);

  this->nlfi_ptr_ = set_nlfi_ptr(dt, u_vect);

  this->RHS->AddDomainIntegrator(this->nlfi_ptr_);
}

/**
 * @brief Configure the Newton solver
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param oper
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void OperatorBase<T, DIM, NLFI, LHS_NLFI>::SetNewtonAlgorithm(mfem::Operator* oper) {
  this->rhs_solver_ =
      new NLSolver(this->nl_solver_, this->nl_solver_params_, this->solver_, this->solver_params_,
                   this->precond_, this->precond_params_, *oper);

  this->newton_solver_ = this->rhs_solver_->get_nl_solver();
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief Compute L2 error
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @tparam LHS_NLFI
 * @param it
 * @param t
 * @param dt
 * @param id_var
 * @param name
 * @param u
 * @param solution_func
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void OperatorBase<T, DIM, NLFI, LHS_NLFI>::ComputeError(
    const int& it, const double& t, const double& dt, const int id_var, const std::string& name,
    const mfem::Vector& u, std::function<double(const mfem::Vector&, double)> solution_func) {
  Catch_Time_Section("OperatorBase::ComputeError");

  mfem::ParGridFunction gf(this->fes_[id_var]);
  mfem::ParGridFunction zero(this->fes_[id_var]);
  zero = 0.0;

  gf.SetFromTrueDofs(u);
  mfem::FunctionCoefficient solution_coef(solution_func);
  solution_coef.SetTime(t);

  // // IR
  // const mfem::FiniteElement *fe = this->fes_[id_var]->GetFE(0);
  // const mfem::IntegrationRule *ir[mfem::Geometry::NumGeom];

  // const mfem::ElementTransformation *Tr = this->fes_[id_var]->GetElementTransformation(0);

  // for (int i = 0; i < mfem::Geometry::NumGeom; ++i) {
  //   std::cout << " i " << i << " Tr->OrderW() " << Tr->OrderW() << " order " << fe->GetOrder()
  //             << std::endl;
  //   ir[i] = &(mfem::IntRules.Get(i, 2 * fe->GetOrder() + 3));
  // }

  const auto errorL2 = gf.ComputeLpError(2., solution_coef);
  const auto errorLinf = gf.ComputeLpError(mfem::infinity(), solution_coef);
  const auto norm_solution = zero.ComputeLpError(2, solution_coef);
  const auto normalized_error = errorL2 / norm_solution;

  this->time_specialized_.emplace(IterationKey(it, dt, t),
                                  SpecializedValue(name + "_L2-error[-]", errorL2));
  this->time_specialized_.emplace(
      IterationKey(it, dt, t),
      SpecializedValue(name + "_L2-error normalized[-]", normalized_error));
  this->time_specialized_.emplace(IterationKey(it, dt, t),
                                  SpecializedValue(name + "_Linf-error [-]", errorLinf));
}

/**
 * @brief Compute integrals of a variable over the domain
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @tparam LHS_NLFI
 * @param it
 * @param t
 * @param dt
 * @param id_var
 * @param name
 * @param u
 * @param integral_threshold
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void OperatorBase<T, DIM, NLFI, LHS_NLFI>::ComputeIntegral(
    const int& it, const double& t, const double& dt, const int id_var, const std::string& name,
    const mfem::Vector& u, const double lower_bound, const double upper_bound) {
  Catch_Time_Section("OperatorBase::ComputeIntegral");

  mfem::ParGridFunction gf(this->fes_[id_var]);
  mfem::Vector u_cut(u.Size());
  std::transform(u.begin(), u.end(), u_cut.begin(), [&](auto value) {
    if (value > upper_bound || value < lower_bound) {
      return 0.0;
    }
    return value;
  });
  gf.SetFromTrueDofs(u_cut);

  double integral = 0.0;
  double domain_volume = 0.0;
  mfem::Vector vals;
  const mfem::FiniteElement* fe;
  mfem::ElementTransformation* Tr;

  for (int i = 0; i < this->fes_[id_var]->GetNE(); i++) {
    fe = this->fes_[id_var]->GetFE(i);
    const mfem::IntegrationRule* ir;

    int intorder = 2 * fe->GetOrder() + 3;  // <----------
    ir = &(mfem::IntRules.Get(fe->GetGeomType(), intorder));

    mfem::real_t int_elem = 0.0;
    mfem::real_t domain_volume_elem = 0.0;

    gf.GetValues(i, *ir, vals);
    Tr = this->fes_[id_var]->GetElementTransformation(i);
    for (int j = 0; j < ir->GetNPoints(); j++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(j);
      Tr->SetIntPoint(&ip);

      int_elem += vals(j) * ip.weight * Tr->Weight();
      domain_volume_elem += ip.weight * Tr->Weight();
    }

    integral += int_elem;
    domain_volume += domain_volume_elem;
  }

  mfem::real_t global_integral = 0.0;
  mfem::real_t global_domain_volume = 0.0;
  MPI_Allreduce(&integral, &global_integral, 1, mfem::MPITypeMap<mfem::real_t>::mpi_type, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(&domain_volume, &global_domain_volume, 1, mfem::MPITypeMap<mfem::real_t>::mpi_type,
                MPI_SUM, MPI_COMM_WORLD);

  integral = global_integral;
  domain_volume = global_domain_volume;

  const double average = integral / domain_volume;
  this->time_specialized_.emplace(IterationKey(it, dt, t),
                                  SpecializedValue(name + "_integral[-]", integral));
  this->time_specialized_.emplace(IterationKey(it, dt, t),
                                  SpecializedValue(name + "_average[-]", average));
}

/**
 * @brief Compute the position of an isovalue
 * @param it current iteration
 * @param t current time
 * @param dt current timestep
 * @param u unknown vector
 * @param iso_value value of the solution
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void OperatorBase<T, DIM, NLFI, LHS_NLFI>::ComputeIsoVal(const int& it, const double& t,
                                                         const double& dt, const int var_id,
                                                         const std::string& var_name,
                                                         const mfem::Vector& u,
                                                         const double& iso_value) {
  Catch_Time_Section("OperatorBase::ComputeIsoVal");
  std::vector<std::string> vstr = {"x", "y", "z"};
  std::multimap<IterationKey, SpecializedValue> variable_time_iso_specialized;

  if (!this->time_iso_specialized_[var_name].empty()) {
    variable_time_iso_specialized.merge(this->time_iso_specialized_[var_name]);
  }

  mfem::ParGridFunction gf(this->fes_[var_id]);
  gf.SetFromTrueDofs(u);
  std::vector<mfem::Vector> iso_points;

  for (int i = 0; i < this->fes_[var_id]->GetNE(); i++) {
    const mfem::FiniteElement* el = this->fes_[var_id]->GetFE(i);
    size_t dim = el->GetDim();
    mfem::Array<int> dofs;
    this->fes_[var_id]->GetElementDofs(i, dofs);

    mfem::DenseMatrix dof_coords;
    mfem::ElementTransformation* Tr = this->fes_[var_id]->GetElementTransformation(i);
    Tr->Transform(el->GetNodes(), dof_coords);
    std::vector<double> dof_val;
    dof_val.reserve(dofs.Size());
    for (int j = 0; j < dofs.Size(); j++) {
      dof_val.emplace_back(u(dofs[j]) - iso_value);
    }

    // At least one positive and one negative value
    bool has_positive_value =
        std::any_of(dof_val.begin(), dof_val.end(), [](double x) { return x > 0; });
    bool has_negative_value =
        std::any_of(dof_val.begin(), dof_val.end(), [](double x) { return x < 0; });

    if (has_positive_value && has_negative_value) {
      for (int j = 0; j < dofs.Size(); j++) {
        mfem::Vector coord1(DIM);
        double val1 = dof_val[j];
        dof_coords.GetColumn(j, coord1);
        for (int k = j + 1; k < dofs.Size(); k++) {
          double val2 = dof_val[k];

          if (val1 * val2 < 0) {
            const double abs_val1 = std::abs(val1);
            const double abs_val2 = std::abs(val2);
            double t = abs_val1 / (abs_val1 + abs_val2);

            mfem::Vector coord2(DIM), iso_coord(DIM);
            iso_coord = mfem::infinity();
            dof_coords.GetColumn(k, coord2);

            for (int d = 0; d < coord1.Size(); d++) {
              iso_coord[d] = coord1[d] + t * (coord2[d] - coord1[d]);
            }
            iso_points.emplace_back(std::move(iso_coord));
          }
        }
      }
      for (size_t i = 0; i < iso_points.size(); i++) {
        for (size_t d = 0; d < DIM; d++) {
          variable_time_iso_specialized.emplace(IterationKey(it, dt, t),
                                                SpecializedValue(vstr[d], iso_points[i](d)));
        }
      }
    }
  }
  this->time_iso_specialized_[var_name] = variable_time_iso_specialized;
}

/**
 * @brief Get the source term by equation of the PDEs
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void OperatorBase<T, DIM, NLFI, LHS_NLFI>::get_source_term(
    const int id_block, const std::function<double(const mfem::Vector&, double)>& src_func,
    mfem::Vector& source_term, mfem::ParLinearForm* RHSS) const {
  mfem::FunctionCoefficient src(src_func);
  src.SetTime(this->current_time_ + this->current_dt_);

  RHSS->AddDomainIntegrator(new mfem::DomainLFIntegrator(src));
  RHSS->Assemble();

  source_term.SetSize(this->fes_[id_block]->GetTrueVSize());
  RHSS->ParallelAssemble(source_term);

  source_term.SetSubVector(this->ess_tdof_list_[id_block], 0.);
}

/**
 * @brief Get a of a multimap of integral values calculated at a given iteration (see
 * computeEnergies and computeError)
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @return const std::map<std::tuple<int, double, double>, double>
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
const std::map<std::string, std::multimap<IterationKey, SpecializedValue>>
OperatorBase<T, DIM, NLFI, LHS_NLFI>::get_time_iso_specialized() const {
  return this->time_iso_specialized_;
}

/**
 * @brief Clear time_iso_specialized_ container
 *
 * @tparam T
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void OperatorBase<T, DIM, NLFI, LHS_NLFI>::clear_iso_time_specialized() {
  this->time_iso_specialized_.clear();
}

/**
 * @brief Get a of a multimap of integral values calculated at a given iteration (see
 * computeEnergies and computeError)
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @return const std::map<std::tuple<int, double, double>, double>
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
const std::multimap<IterationKey, SpecializedValue>
OperatorBase<T, DIM, NLFI, LHS_NLFI>::get_time_specialized() const {
  return this->time_specialized_;
}

/**
 * @brief Clear time_specialized_ container
 *
 * @tparam T
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void OperatorBase<T, DIM, NLFI, LHS_NLFI>::clear_time_specialized() {
  this->time_specialized_.clear();
}

/**
 * @brief Set the default options for the nonlinear algorithm and associated solvers
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void OperatorBase<T, DIM, NLFI, LHS_NLFI>::set_default_solver() {
  auto nl_params =
      Parameters(Parameter("description", "Newton Algorithm"), Parameter("iterative_mode", false));
  auto s_params = Parameters(Parameter("description", "Default Solver for Newton Algorithm"));
  auto p_params =
      Parameters(Parameter("description", "Default Preconditioner for Newton Algorithm"));

  this->nl_solver_ = NLSolverType::NEWTON;
  this->nl_solver_params_ = nl_params;
  this->solver_ = HypreSolverType::HYPRE_GMRES;
  this->solver_params_ = s_params;
  this->precond_ = HyprePreconditionerType::HYPRE_ILU;
  this->precond_params_ = p_params;
}

/**
 * @brief Overload the nonlinear algorithm
 * @remark NewtonSolver is the only one implemented
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param NLSOLVER
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void OperatorBase<T, DIM, NLFI, LHS_NLFI>::overload_nl_solver(NLSolverType NLSOLVER) {
  this->nl_solver_ = NLSOLVER;
}

/**
 * @brief Overload Parameters associated with nonlinear algorithm
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param NLSOLVER
 * @param nl_params
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void OperatorBase<T, DIM, NLFI, LHS_NLFI>::overload_nl_solver(NLSolverType NLSOLVER,
                                                              const Parameters& nl_params) {
  this->nl_solver_ = NLSOLVER;
  this->nl_solver_params_ = nl_params;
}

/**
 * @brief  Overload the default linear solver used for the LHS
 * @remark only used for explicit time-scheme
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param SOLVER
 * @param s_params
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void OperatorBase<T, DIM, NLFI, LHS_NLFI>::overload_solver(VSolverType SOLVER) {
  this->solver_ = SOLVER;
}

/**
 * @brief  Overload the parameter for the linear solver used for the LHS
 * @remark only used for explicit time-scheme
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param SOLVER
 * @param s_params
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void OperatorBase<T, DIM, NLFI, LHS_NLFI>::overload_solver(VSolverType SOLVER,
                                                           const Parameters& s_params) {
  this->solver_ = SOLVER;
  this->solver_params_ = s_params;
}

/**
 * @brief Overload the preconditioner used by the solver in the NL algorithm.
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param PRECOND
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void OperatorBase<T, DIM, NLFI, LHS_NLFI>::overload_preconditioner(VSolverType PRECOND) {
  this->precond_ = PRECOND;
}

/**
 * @brief  Overload the parameters of the preconditioner used by the solver in the NL algorithm.
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param PRECOND
 * @param p_params
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void OperatorBase<T, DIM, NLFI, LHS_NLFI>::overload_preconditioner(VSolverType PRECOND,
                                                                   const Parameters& p_params) {
  this->precond_ = PRECOND;
  this->precond_params_ = p_params;
}

/**
 * @brief Destroy the Operator Base<T, DIM,  NLFI>:: Operator Base object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
OperatorBase<T, DIM, NLFI, LHS_NLFI>::~OperatorBase() {
  if (this->nlfi_ptr_ != nullptr) {
    delete this->nlfi_ptr_;
  }
}

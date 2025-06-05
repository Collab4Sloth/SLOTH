/**
 * @file OperatorBase.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Base class for building Steady and TimeDependent PhaseFieldOperators
 * @version 0.1
 * @date 2024-07-25
 *
 * Copyright CEA (c) 2024
 *
 */
#include <algorithm>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "AnalyticalFunctions/AnalyticalFunctions.hpp"
#include "BCs/BoundaryConditions.hpp"
#include "Coefficients/EnergyCoefficient.hpp"
#include "Coefficients/MobilityCoefficient.hpp"
#include "Coefficients/PhaseChangeFunction.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Profiling/Profiling.hpp"
#include "Solvers/LSolver.hpp"
#include "Solvers/NLSolver.hpp"
#include "Spatial/Spatial.hpp"
#include "Utils/Utils.hpp"
#include "Variables/Variable.hpp"
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
template <class T, int DIM, class NLFI>
class OperatorBase : public mfem::Operator {
 private:
  T *fecollection_;

  NLSolverType nl_solver_;
  VSolverType solver_;
  VSolverType precond_;
  Parameters nl_solver_params_;
  Parameters solver_params_;
  Parameters precond_params_;
  void set_default_solver();

 protected:
  std::vector<Variables<T, DIM> *> auxvariables_;
  std::string description_{"UNKNOWN OPERATOR"};
  const Parameter default_p_ = Parameter("default parameter", false);
  const Parameters default_params_ = Parameters(default_p_);
  const Parameters &params_;
  /// Time integral results
  std::multimap<IterationKey, SpecializedValue> time_specialized_, time_iso_specialized_;
  /// Finite Element Space and BCs
  mfem::Array<mfem::ParFiniteElementSpace *> fes_;
  mfem::ParFiniteElementSpace *fespace_;
  std::vector<mfem::Array<int>> ess_tdof_list_;
  mfem::Array<int> block_trueOffsets_;

  /// Left-Hand-Side + Right-Hand-Side
  mfem::ParBlockNonlinearForm *LHS;
  mfem::ParBlockNonlinearForm *N;

  NLSolver *rhs_solver_;
  std::shared_ptr<mfem::NewtonSolver> newton_solver_;

  /// Boundary conditions
  std::vector<BoundaryConditions<T, DIM> *> bcs_;

  std::function<double(const mfem::Vector &, double)> src_func_;

  double current_dt_;
  double current_time_;
  int height_;
  mutable mfem::Vector z;  // auxiliary vector
  NLFI *nlfi_ptr_;

  void build_lhs_nonlinear_form(const double dt, const std::vector<mfem::Vector> &u);
  void build_nonlinear_form(const double dt, const std::vector<mfem::Vector> &u);
  void SetNewtonAlgorithm(mfem::Operator *oper);

  int compute_total_size(const std::vector<SpatialDiscretization<T, DIM> *> &spatials);

 public:
  explicit OperatorBase(std::vector<SpatialDiscretization<T, DIM> *> spatials);

  OperatorBase(std::vector<SpatialDiscretization<T, DIM> *> spatials,
               AnalyticalFunctions<DIM> source_term_name);
  OperatorBase(std::vector<SpatialDiscretization<T, DIM> *> spatials, const Parameters &params);

  OperatorBase(std::vector<SpatialDiscretization<T, DIM> *> spatials, const Parameters &params,
               AnalyticalFunctions<DIM> source_term_name);

  void ComputeError(const int &it, const double &t, const double &dt, const mfem::Vector &u,
                    std::function<double(const mfem::Vector &, double)> solution_func);
  void ComputeIsoVal(const int &it, const double &t, const double &dt, const mfem::Vector &u,
                     const double &iso_value);
  void get_source_term(mfem::Vector &source_term, mfem::ParLinearForm *RHHS) const;

  const std::multimap<IterationKey, SpecializedValue> get_time_specialized() const;
  const std::multimap<IterationKey, SpecializedValue> get_time_iso_specialized() const;

  virtual ~OperatorBase();

  std::string get_description() { return this->description_; }

  // User-defined Solvers
  void overload_nl_solver(NLSolverType NLSOLVER);
  void overload_nl_solver(NLSolverType NLSOLVER, const Parameters &nl_params);

  void overload_solver(VSolverType SOLVER);
  void overload_solver(VSolverType SOLVER, const Parameters &s_params);

  void overload_preconditioner(VSolverType PRECOND);
  void overload_preconditioner(VSolverType PRECOND, const Parameters &p_params);

  // Virtual methods
  // virtual void initialize(const double &initial_time, Variables<T, DIM> &vars);
  virtual void initialize(const double &initial_time, Variables<T, DIM> &vars,
                          std::vector<Variables<T, DIM> *> auxvars);

  // Pure virtual methods
  virtual void set_default_properties() = 0;
  virtual void SetConstantParameters(const double dt, const std::vector<mfem::Vector> &u_vect) = 0;
  virtual void SetTransientParameters(const double dt, const std::vector<mfem::Vector> &u_vect) = 0;
  virtual void solve(std::vector<std::unique_ptr<mfem::Vector>> &vect_unk, double &next_time,
                     const double &current_time, double current_time_step, const int iter) = 0;
  virtual NLFI *set_nlfi_ptr(const double dt, const std::vector<mfem::Vector> &u) = 0;
  virtual void get_parameters() = 0;
  virtual void ComputeEnergies(const int &it, const double &dt, const double &t,
                               const mfem::Vector &u) = 0;
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

template <class T, int DIM, class NLFI>
int OperatorBase<T, DIM, NLFI>::compute_total_size(
    const std::vector<SpatialDiscretization<T, DIM> *> &spatials) {
  int total_size = 0;
  for (const auto *s : spatials) {
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
template <class T, int DIM, class NLFI>
OperatorBase<T, DIM, NLFI>::OperatorBase(std::vector<SpatialDiscretization<T, DIM> *> spatials)
    : mfem::Operator(this->compute_total_size(spatials)),
      params_(default_params_),
      N(NULL),
      LHS(NULL),
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
  for (const auto *s : spatials) {
    this->fes_[i] = s->get_finite_element_space();
    this->block_trueOffsets_[i + 1] = this->fes_[i]->GetTrueVSize();
    i++;
  }
  this->block_trueOffsets_.PartialSum();
  this->fespace_ = this->fes_[0];
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
template <class T, int DIM, class NLFI>
OperatorBase<T, DIM, NLFI>::OperatorBase(std::vector<SpatialDiscretization<T, DIM> *> spatials,
                                         AnalyticalFunctions<DIM> source_term_name)
    : mfem::Operator(this->compute_total_size(spatials)),
      params_(default_params_),
      N(NULL),
      LHS(NULL),
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
  for (const auto *s : spatials) {
    this->fes_[i] = s->get_finite_element_space();
    this->block_trueOffsets_[i + 1] = this->fes_[i]->GetTrueVSize();
    i++;
  }
  this->block_trueOffsets_.PartialSum();
  this->fespace_ = this->fes_[0];
  this->src_func_ = source_term_name.getFunction();
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
template <class T, int DIM, class NLFI>
OperatorBase<T, DIM, NLFI>::OperatorBase(std::vector<SpatialDiscretization<T, DIM> *> spatials,
                                         const Parameters &params)
    : mfem::Operator(this->compute_total_size(spatials)),
      params_(params),
      N(NULL),
      LHS(NULL),
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
  for (const auto *s : spatials) {
    this->fes_[i] = s->get_finite_element_space();
    this->block_trueOffsets_[i + 1] = this->fes_[i]->GetTrueVSize();
    i++;
  }
  this->block_trueOffsets_.PartialSum();
  this->fespace_ = this->fes_[0];
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
template <class T, int DIM, class NLFI>
OperatorBase<T, DIM, NLFI>::OperatorBase(std::vector<SpatialDiscretization<T, DIM> *> spatials,
                                         const Parameters &params,
                                         AnalyticalFunctions<DIM> source_term_name)
    : mfem::Operator(this->compute_total_size(spatials)),
      params_(params),
      N(NULL),
      LHS(NULL),
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
  for (const auto *s : spatials) {
    this->fes_[i] = s->get_finite_element_space();
    this->block_trueOffsets_[i + 1] = this->fes_[i]->GetTrueVSize();
    i++;
  }
  this->block_trueOffsets_.PartialSum();
  this->fespace_ = this->fes_[0];
  this->src_func_ = source_term_name.getFunction();
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
template <class T, int DIM, class NLFI>
void OperatorBase<T, DIM, NLFI>::initialize(const double &initial_time, Variables<T, DIM> &vars,
                                            std::vector<Variables<T, DIM> *> auxvars) {
  Catch_Time_Section("OperatorBase::initialize");

  this->auxvariables_ = auxvars;
  const auto nvars = vars.get_variables_number();
  std::vector<mfem::Vector> u_vect;
  u_vect.reserve(nvars);
  for (auto iv = 0; iv < nvars; iv++) {
    auto &vv = vars.getIVariable(iv);
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
 * @brief Build the NonLinear Form Integrator
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param dt
 * @param u
 */
template <class T, int DIM, class NLFI>
void OperatorBase<T, DIM, NLFI>::build_nonlinear_form(const double dt,
                                                      const std::vector<mfem::Vector> &u_vect) {
  ////////////////////////////////////////////
  // PhaseField non linear form
  ////////////////////////////////////////////

  // TODO(cci) change methods

  auto u = u_vect[0];

  if (N != nullptr) {
    delete N;
  }
  N = new mfem::ParBlockNonlinearForm(this->fes_);
  // mfem::ParGridFunction un_gf(this->fespace_);
  // un_gf.SetFromTrueDofs(u);
  // mfem::ParGridFunction un(this->fespace_);
  // un.SetFromTrueDofs(u);

  // TODO(cci) : pass a vector of GF.
  this->nlfi_ptr_ = set_nlfi_ptr(dt, u_vect);

  N->AddDomainIntegrator(this->nlfi_ptr_);
  // TODO(cci) check BCs if always usefull here
  // N->SetEssentialTrueDofs(this->ess_tdof_list_[0]);
}

/**
 * @brief Build the LHS
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param dt
 * @param u
 */
template <class T, int DIM, class NLFI>
void OperatorBase<T, DIM, NLFI>::build_lhs_nonlinear_form(const double dt,
                                                          const std::vector<mfem::Vector> &u) {
  ////////////////////////////////////////////
  // PhaseField non linear form : LHS
  ////////////////////////////////////////////

  // auto u = u_vect[0];

  if (LHS != nullptr) {
    delete LHS;
  }
  LHS = new mfem::ParBlockNonlinearForm(this->fes_);
  // mfem::ParGridFunction un_gf(this->fespace_);
  // un_gf.SetFromTrueDofs(u);
  // mfem::ParGridFunction un(this->fespace_);
  // un.SetFromTrueDofs(u);

  // TODO(cci) : pass a vector of GF.

  const Parameters &all_params = this->params_ - this->default_p_;
  std::vector<mfem::ParGridFunction> vun;
  for (int i = 0; i < u.size(); i++) {
    mfem::ParGridFunction un(this->fes_[i]);
    un.SetFromTrueDofs(u[i]);
    vun.emplace_back(un);
  }
  TimeNLFormIntegrator<Variables<T, DIM>> *time_nlfi_ptr =
      new TimeNLFormIntegrator<Variables<T, DIM>>(vun, all_params, this->auxvariables_);

  LHS->AddDomainIntegrator(time_nlfi_ptr);

  // TODO(cci) check BCs
  // N->SetEssentialTrueDofs(this->ess_tdof_list_[0]);
}

/**
 * @brief Configure the Newton solver
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param oper
 */
template <class T, int DIM, class NLFI>
void OperatorBase<T, DIM, NLFI>::SetNewtonAlgorithm(mfem::Operator *oper) {
  ////////////////////////////////////////////
  // Newton Solver
  ////////////////////////////////////////////
  // if (this->precond_ != PreconditionerType::NO) {
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
 * @param u unknown vector
 * @return const double
 */
template <class T, int DIM, class NLFI>
void OperatorBase<T, DIM, NLFI>::ComputeError(
    const int &it, const double &t, const double &dt, const mfem::Vector &u,
    std::function<double(const mfem::Vector &, double)> solution_func) {
  Catch_Time_Section("OperatorBase::ComputeError");

  mfem::ParGridFunction gf(this->fespace_);
  mfem::ParGridFunction zero(this->fespace_);
  zero = 0.0;

  gf.SetFromTrueDofs(u);
  mfem::FunctionCoefficient solution_coef(solution_func);
  solution_coef.SetTime(t);

  const auto errorL2 = gf.ComputeLpError(2., solution_coef);
  const auto errorLinf = gf.ComputeLpError(mfem::infinity(), solution_coef);
  const auto norm_solution = zero.ComputeLpError(2, solution_coef);
  const auto normalized_error = errorL2 / norm_solution;

  this->time_specialized_.emplace(IterationKey(it, dt, t),
                                  SpecializedValue("L2-error[-]", errorL2));
  this->time_specialized_.emplace(IterationKey(it, dt, t),
                                  SpecializedValue("L2-error normalized[-]", normalized_error));
  this->time_specialized_.emplace(IterationKey(it, dt, t),
                                  SpecializedValue("Linf-error [-]", errorLinf));
}

/**
 * @brief Compute the position of an isovalue
 * @param it current iteration
 * @param t current time
 * @param dt current timestep
 * @param u unknown vector
 * @param iso_value value of the solution
 */
template <class T, int DIM, class NLFI>
void OperatorBase<T, DIM, NLFI>::ComputeIsoVal(const int &it, const double &t, const double &dt,
                                               const mfem::Vector &u, const double &iso_value) {
  Catch_Time_Section("OperatorBase::ComputeIsoVal");
  std::vector<std::string> vstr = {"x", "y", "z"};

  mfem::ParGridFunction gf(this->fespace_);
  gf.SetFromTrueDofs(u);
  std::vector<mfem::Vector> iso_points;

  for (int i = 0; i < this->fespace_->GetNE(); i++) {
    const mfem::FiniteElement *el = this->fespace_->GetFE(i);
    size_t dim = el->GetDim();
    mfem::Array<int> dofs;
    this->fespace_->GetElementDofs(i, dofs);

    mfem::DenseMatrix dof_coords;
    mfem::ElementTransformation *Tr = this->fespace_->GetElementTransformation(i);
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
          this->time_iso_specialized_.emplace(IterationKey(it, dt, t),
                                              SpecializedValue(vstr[d], iso_points[i](d)));
        }
      }
    }
  }
}

/**
 * @brief get the source term
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 */
template <class T, int DIM, class NLFI>
void OperatorBase<T, DIM, NLFI>::get_source_term(mfem::Vector &source_term,
                                                 mfem::ParLinearForm *RHSS) const {
  mfem::FunctionCoefficient src(this->src_func_);

  src.SetTime(this->current_time_);
  RHSS->AddDomainIntegrator(new mfem::DomainLFIntegrator(src));
  RHSS->Assemble();

  // BCs
  source_term.SetSize(this->fespace_->GetTrueVSize());
  RHSS->ParallelAssemble(source_term);

  source_term.SetSubVector(this->ess_tdof_list_[0], 0.);
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
template <class T, int DIM, class NLFI>
const std::multimap<IterationKey, SpecializedValue>
OperatorBase<T, DIM, NLFI>::get_time_iso_specialized() const {
  return this->time_iso_specialized_;
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
template <class T, int DIM, class NLFI>
const std::multimap<IterationKey, SpecializedValue>
OperatorBase<T, DIM, NLFI>::get_time_specialized() const {
  return this->time_specialized_;
}

/**
 * @brief Set the default options for the non linear algorithm and associated solvers
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 */
template <class T, int DIM, class NLFI>
void OperatorBase<T, DIM, NLFI>::set_default_solver() {
  auto nl_params =
      Parameters(Parameter("description", "Newton Algorithm"), Parameter("iterative_mode", false));
  auto s_params = Parameters(Parameter("description", "HypreGMRES Solver"));
  auto p_params = Parameters(Parameter("description", "HypreILU preconditioner"));

  this->nl_solver_ = NLSolverType::NEWTON;
  this->nl_solver_params_ = nl_params;
  this->solver_ = HypreSolverType::HYPRE_GMRES;
  this->solver_params_ = s_params;
  this->precond_ = HyprePreconditionerType::HYPRE_ILU;
  this->precond_params_ = p_params;
}

/**
 * @brief Overload  the non linear algorithm
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param NLSOLVER
 */
template <class T, int DIM, class NLFI>
void OperatorBase<T, DIM, NLFI>::overload_nl_solver(NLSolverType NLSOLVER) {
  this->nl_solver_ = NLSOLVER;
}

/**
 * @brief Overload  the non linear algorithm with Parameters
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param NLSOLVER
 * @param nl_params
 */
template <class T, int DIM, class NLFI>
void OperatorBase<T, DIM, NLFI>::overload_nl_solver(NLSolverType NLSOLVER,
                                                    const Parameters &nl_params) {
  this->nl_solver_ = NLSOLVER;

  this->nl_solver_params_ = nl_params;
}

/**
 * @brief  Overload the default linear solver
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param SOLVER
 * @param s_params
 */
template <class T, int DIM, class NLFI>
void OperatorBase<T, DIM, NLFI>::overload_solver(VSolverType SOLVER) {
  this->solver_ = SOLVER;
}

/**
 * @brief   Overload the default linear solver with parameters
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param SOLVER
 * @param s_params
 */
template <class T, int DIM, class NLFI>
void OperatorBase<T, DIM, NLFI>::overload_solver(VSolverType SOLVER, const Parameters &s_params) {
  this->solver_ = SOLVER;
  this->solver_params_ = s_params;
}

/**
 * @brief Overload the preconditioner
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param PRECOND
 */
template <class T, int DIM, class NLFI>
void OperatorBase<T, DIM, NLFI>::overload_preconditioner(VSolverType PRECOND) {
  this->precond_ = PRECOND;
}

/**
 * @brief Overload the preconditioner with parameters
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param PRECOND
 * @param p_params
 */
template <class T, int DIM, class NLFI>
void OperatorBase<T, DIM, NLFI>::overload_preconditioner(VSolverType PRECOND,
                                                         const Parameters &p_params) {
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
template <class T, int DIM, class NLFI>
OperatorBase<T, DIM, NLFI>::~OperatorBase() {
  if (this->nlfi_ptr_ != nullptr) {
    delete this->nlfi_ptr_;
  }
}

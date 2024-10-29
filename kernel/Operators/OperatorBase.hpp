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
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "BCs/BoundaryConditions.hpp"
#include "Coefficients/EnergyCoefficient.hpp"
#include "Coefficients/MobilityCoefficient.hpp"
#include "Coefficients/PhaseChangeFunction.hpp"
#include "Coefficients/SourceTermCoefficient.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Profiling/Profiling.hpp"
#include "Solvers/LSolver.hpp"
#include "Solvers/NLSolver.hpp"
#include "Spatial/Spatial.hpp"
#include "Utils/AnalyticalFunctions.hpp"
#include "Utils/PhaseFieldConstants.hpp"
#include "Utils/PhaseFieldOptions.hpp"
#include "Utils/UtilsForDebug.hpp"
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

  const Parameters default_params_ = Parameters(Parameter("default parameter", false));
  const Parameters &params_;
  /// Time integral results
  std::multimap<IterationKey, SpecializedValue> time_specialized_;

  /// Finite Element Space and BCs
  mfem::ParFiniteElementSpace *fespace_;
  mfem::Array<int> ess_tdof_list_;

  /// Right-Hand-Side
  mfem::ParNonlinearForm *N;
  NLSolver *rhs_solver_;
  std::shared_ptr<mfem::NewtonSolver> newton_solver_;

  /// Boundary conditions
  BoundaryConditions<T, DIM> *bcs_;

  std::function<double(const mfem::Vector &, double)> src_func_;

  double current_dt_;
  double current_time_;
  int height_;
  mutable mfem::Vector z;  // auxiliary vector
  NLFI *nlfi_ptr_;
  void build_nonlinear_form(const double dt, const mfem::Vector &u);
  void SetNewtonAlgorithm(mfem::Operator *oper);

  std::vector<mfem::ParGridFunction> get_auxiliary_gf();

 public:
  explicit OperatorBase(SpatialDiscretization<T, DIM> const *spatial);

  OperatorBase(SpatialDiscretization<T, DIM> const *spatial,
               AnalyticalFunctions<DIM> source_term_name);
  OperatorBase(SpatialDiscretization<T, DIM> const *spatial, const Parameters &params);

  OperatorBase(SpatialDiscretization<T, DIM> const *spatial, const Parameters &params,
               AnalyticalFunctions<DIM> source_term_name);

  void ComputeError(const int &it, const double &t, const double &dt, const mfem::Vector &u,
                    std::function<double(const mfem::Vector &, double)> solution_func);
  void get_source_term(mfem::Vector &source_term, mfem::ParLinearForm *RHHS) const;

  const std::multimap<IterationKey, SpecializedValue> get_time_specialized() const;

  virtual ~OperatorBase();

  std::string get_description() { return this->description_; }

  // User-defined Solvers
  void overload_nl_solver(NLSolverType NLSOLVER, const Parameters &nl_params);
  void overload_solver(VSolverType SOLVER, const Parameters &s_params);
  void overload_solver(VSolverType SOLVER, const Parameters &s_params, VSolverType PRECOND,
                       const Parameters &p_params);
  void overload_preconditioner(VSolverType PRECOND, const Parameters &p_params);

  // Virtual methods
  // virtual void initialize(const double &initial_time, Variables<T, DIM> &vars);
  virtual void initialize(const double &initial_time, Variables<T, DIM> &vars,
                          std::vector<Variables<T, DIM> *> auxvars);

  // Pure virtual methods
  virtual void set_default_properties() = 0;
  virtual void SetConstantParameters(const double dt, mfem::Vector &u) = 0;
  virtual void SetTransientParameters(const double dt, const mfem::Vector &u) = 0;
  virtual void solve(mfem::Vector &unk, double &next_time, const double &current_time,
                     double current_time_step, const int iter) = 0;
  virtual NLFI *set_nlfi_ptr(const double dt, const mfem::Vector &u) = 0;
  virtual void get_parameters() = 0;
  virtual void ComputeEnergies(const int &it, const double &dt, const double &t,
                               const mfem::Vector &u) = 0;
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Construct a new OperatorBase object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param spatial
 */
template <class T, int DIM, class NLFI>
OperatorBase<T, DIM, NLFI>::OperatorBase(SpatialDiscretization<T, DIM> const *spatial)
    : mfem::Operator(spatial->getSize()),
      params_(default_params_),
      N(NULL),
      current_dt_(0.0),
      current_time_(0.0),
      height_(height),
      z(height),
      nlfi_ptr_(nullptr) {
  this->fespace_ = spatial->get_finite_element_space();
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
OperatorBase<T, DIM, NLFI>::OperatorBase(SpatialDiscretization<T, DIM> const *spatial,
                                         AnalyticalFunctions<DIM> source_term_name)
    : mfem::Operator(spatial->getSize()),
      params_(default_params_),
      N(NULL),
      current_dt_(0.0),
      current_time_(0.0),
      height_(height),
      z(height),
      nlfi_ptr_(nullptr) {
  this->fespace_ = spatial->get_finite_element_space();
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
OperatorBase<T, DIM, NLFI>::OperatorBase(SpatialDiscretization<T, DIM> const *spatial,
                                         const Parameters &params)
    : mfem::Operator(spatial->getSize()),
      params_(params),
      N(NULL),
      current_dt_(0.0),
      current_time_(0.0),
      height_(height),
      z(height),
      nlfi_ptr_(nullptr) {
  this->fespace_ = spatial->get_finite_element_space();
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
OperatorBase<T, DIM, NLFI>::OperatorBase(SpatialDiscretization<T, DIM> const *spatial,
                                         const Parameters &params,
                                         AnalyticalFunctions<DIM> source_term_name)
    : mfem::Operator(spatial->getSize()),
      params_(params),
      N(NULL),
      current_dt_(0.0),
      current_time_(0.0),
      height_(height),
      z(height),
      nlfi_ptr_(nullptr) {
  this->fespace_ = spatial->get_finite_element_space();
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

  auto &vv = vars.getIVariable(0);
  auto u = vv.get_unknown();

  this->bcs_ = vv.get_boundary_conditions();
  this->ess_tdof_list_ = this->bcs_->GetEssentialDofs();
  this->bcs_->SetBoundaryConditions(u);

  this->SetConstantParameters(this->current_dt_, u);

  this->SetTransientParameters(this->current_dt_, u);

  vv.update(u);
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
void OperatorBase<T, DIM, NLFI>::build_nonlinear_form(const double dt, const mfem::Vector &u) {
  ////////////////////////////////////////////
  // PhaseField non linear form
  ////////////////////////////////////////////
  if (N != nullptr) {
    delete N;
  }
  N = new mfem::ParNonlinearForm(this->fespace_);
  mfem::ParGridFunction un_gf(this->fespace_);
  un_gf.SetFromTrueDofs(u);
  mfem::ParGridFunction un(this->fespace_);
  un.SetFromTrueDofs(u);

  // if (this->nlfi_ptr_ != nullptr) {
  //   delete this->nlfi_ptr_;
  // }
  // NLFI *this->nlfi_ptr_  = set_this->nlfi_ptr_ (dt, u);
  this->nlfi_ptr_ = set_nlfi_ptr(dt, u);

  N->AddDomainIntegrator(this->nlfi_ptr_);
  N->SetEssentialTrueDofs(this->ess_tdof_list_);
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
  // } else {
  //   this->rhs_solver_ = new NLSolver(this->nl_solver_, this->nl_solver_params_, this->solver_,
  //                                    this->solver_params_, *oper);
  // }
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
  gf.SetFromTrueDofs(u);
  mfem::FunctionCoefficient solution_coef(solution_func);
  solution_coef.SetTime(t);

  const auto error = gf.ComputeLpError(2, solution_coef);

  this->time_specialized_.emplace(IterationKey(it, dt, t), SpecializedValue("L2-error[-]", error));
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

  // source_term = *RHSS.get();
  source_term.SetSubVector(this->ess_tdof_list_, 0.);

  // this->RHS = std::make_unique<mfem::ParLinearForm>(this->fespace_);
  // mfem::FunctionCoefficient src(this->src_func_);

  // src.SetTime(this->current_time_);
  // this->RHS->AddDomainIntegrator(new mfem::DomainLFIntegrator(src));
  // this->RHS->Assemble();

  // // BCs
  // source_term.SetSize(this->fespace_->GetTrueVSize());
  // this->RHS->ParallelAssemble(source_term);

  // source_term.SetSubVector(this->ess_tdof_list_, 0.);
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
 * @brief Overload the default options for the non linear algorithm
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
 * @brief  Overload the default options for solvers
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
 * @brief  Overload the default options for solvers with preconditionners
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param PRECOND
 * @param p_params
 */
template <class T, int DIM, class NLFI>
void OperatorBase<T, DIM, NLFI>::overload_solver(VSolverType SOLVER, const Parameters &s_params,
                                                 VSolverType PRECOND, const Parameters &p_params) {
  this->overload_solver(SOLVER, s_params);
  this->precond_ = PRECOND;
  this->precond_params_ = p_params;
}

/**
 * @brief  Overload the default options for preconditionners
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
 * @brief Return the vector of grid functions associated with the auxiliary variables
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 */
template <class T, int DIM, class NLFI>
std::vector<mfem::ParGridFunction> OperatorBase<T, DIM, NLFI>::get_auxiliary_gf() {
  std::vector<mfem::ParGridFunction> aux_gf;
  if (this->auxvariables_.size() > 0) {
    for (const auto &auxvar_vec : this->auxvariables_) {
      for (const auto &auxvar : auxvar_vec->getVariables()) {
      auto gf = auxvar.get_gf();
      aux_gf.emplace_back(gf);
      }
    }
  }
  return aux_gf;
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

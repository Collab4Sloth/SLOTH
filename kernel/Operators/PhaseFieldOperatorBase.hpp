/*
 * Copyright © CEA 2022
 *
 * \brief PhaseField time-dependent operator built on the basis of ex16.cpp from MFEM repository
 *
 *   After spatial discretization, the phasefield model can be written as:
 *      dphi/dt = M^{-1}(-Kphi)
 *   where phi denotes the phase-field variable and K the phase-field operator (that does not
 *   depends on time)
 *
 * \file PhaseFieldOperatorBase.hpp
 * \author ci230846
 * \date 12/01/2022
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
#include "Operators/ReducedOperator.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Profiling/Profiling.hpp"
#include "Solvers/Solver.hpp"
#include "Spatial/Spatial.hpp"
#include "Utils/AnalyticalFunctions.hpp"
#include "Utils/PhaseFieldConstants.hpp"
#include "Utils/PhaseFieldOptions.hpp"
#include "Utils/UtilsForDebug.hpp"
#include "Variables/Variable.hpp"
#include "Variables/Variables.hpp"
#include "mfem.hpp"

#pragma once

/**
 * @brief PhaseFieldOperatorBase class
 *
 */
template <class T, int DIM, class NLFI>
class PhaseFieldOperatorBase : public mfem::TimeDependentOperator {
 private:
  T *fecollection_;

 protected:
  // Results
  std::multimap<IterationKey, SpecializedValue> time_specialized_;

  mfem::FiniteElementSpace *fespace_;
  mfem::Array<int> ess_tdof_list_;

  /// Mass operator
  mfem::BilinearForm *M;  // mass operator
  Solver *mass_solver_;
  std::shared_ptr<mfem::IterativeSolver>
      M_solver_;  // Krylov solver for inverting the mass matrix M

  mfem::SparseMatrix Mmat;

  /// Right-Hand-Side
  mfem::LinearForm *RHS;
  mfem::NonlinearForm *N;
  Solver *rhs_solver_;
  std::shared_ptr<mfem::IterativeSolver> newton_solver_;

  /// Boundary conditions
  BoundaryConditions<T, DIM> *bcs_;
  Variables<T, DIM> vars_;
  Variables<T, DIM> auxvars_;

  /** Nonlinear operator defining the reduced backward Euler equation for the
      velocity. Used in the implementation of method ImplicitSolve. */
  PhaseFieldReducedOperator *reduced_oper;

  std::function<double(const mfem::Vector &, double)> src_func_;

  double current_dt_;
  mutable mfem::Vector z;  // auxiliary vector
  std::string source_term_name_;

 public:
  PhaseFieldOperatorBase(SpatialDiscretization<T, DIM> *spatial, const Parameters &params,
                         Variables<T, DIM> &vars);
  PhaseFieldOperatorBase(SpatialDiscretization<T, DIM> *spatial, const Parameters &params,
                         Variables<T, DIM> &vars, Variables<T, DIM> &auxvars);

  PhaseFieldOperatorBase(SpatialDiscretization<T, DIM> *spatial, const Parameters &params,
                         Variables<T, DIM> &vars, AnalyticalFunctions<DIM> source_term_name);

  PhaseFieldOperatorBase(SpatialDiscretization<T, DIM> *spatial, const Parameters &params,
                         Variables<T, DIM> &vars, Variables<T, DIM> &auxvars,
                         AnalyticalFunctions<DIM> source_term_name);
  virtual void Mult(const mfem::Vector &u, mfem::Vector &du_dt) const;
  /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
      This is the only requirement for high-order SDIRK implicit integration.*/
  virtual void ImplicitSolve(const double dt, const mfem::Vector &u, mfem::Vector &k);
  void SetConstantParameters(const double dt, mfem::Vector &u);
  void SetTransientParameters(const double dt, const mfem::Vector &u);
  void initialize(const double &initial_time);
  void ComputeError(const int &it, const double &dt, const double &t, const mfem::Vector &u,
                    std::function<double(const mfem::Vector &, double)> solution_func);
  void get_source_term(mfem::Vector &source_term) const;

  const std::multimap<IterationKey, SpecializedValue> get_time_specialized() const;

  virtual ~PhaseFieldOperatorBase();

  // Pure virtual methods
  virtual NLFI *set_nlfi_ptr(const double dt, const mfem::Vector &u) = 0;
  virtual void get_parameters(const Parameters &vectr_param) = 0;
  virtual void ComputeEnergies(const int &it, const double &dt, const double &t,
                               const mfem::Vector &u) = 0;
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Construct a new Phase Field Operator Base< T,  DIM, NLFI>:: Phase Field Operator
 * Base object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param spatial
 * @param params
 * @param vars
 * @param auxvars
 */
template <class T, int DIM, class NLFI>
PhaseFieldOperatorBase<T, DIM, NLFI>::PhaseFieldOperatorBase(SpatialDiscretization<T, DIM> *spatial,
                                                             const Parameters &params,
                                                             Variables<T, DIM> &vars,
                                                             Variables<T, DIM> &auxvars)
    : mfem::TimeDependentOperator(spatial->getSize(), 0.0),
      M(NULL),
      N(NULL),
      reduced_oper(NULL),
      vars_(vars),
      auxvars_(auxvars),
      current_dt_(0.0),
      z(height) {
  this->fespace_ = spatial->get_finite_element_space();
}

/**
 * @brief Construct a new Phase Field Operator Base< T,  DIM, NLFI>:: Phase Field Operator
 * Base object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param spatial
 * @param params
 * @param vars
 */
template <class T, int DIM, class NLFI>
PhaseFieldOperatorBase<T, DIM, NLFI>::PhaseFieldOperatorBase(SpatialDiscretization<T, DIM> *spatial,
                                                             const Parameters &params,
                                                             Variables<T, DIM> &vars)
    : mfem::TimeDependentOperator(spatial->getSize(), 0.0),
      M(NULL),
      N(NULL),
      reduced_oper(NULL),
      vars_(vars),
      current_dt_(0.0),
      z(height) {
  this->fespace_ = spatial->get_finite_element_space();
}
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Construct a new Phase Field Operator Base< T,  DIM, NLFI>:: Phase Field Operator
 * Base object
 *
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
PhaseFieldOperatorBase<T, DIM, NLFI>::PhaseFieldOperatorBase(
    SpatialDiscretization<T, DIM> *spatial, const Parameters &params, Variables<T, DIM> &vars,
    Variables<T, DIM> &auxvars, AnalyticalFunctions<DIM> source_term_name)
    : mfem::TimeDependentOperator(spatial->getSize(), 0.0),
      M(NULL),
      N(NULL),
      reduced_oper(NULL),
      vars_(vars),
      auxvars_(auxvars),
      current_dt_(0.0),
      z(height) {
  this->fespace_ = spatial->get_finite_element_space();
  this->src_func_ = source_term_name.getFunction();

  // auto &vv = vars.get_variable("phi");
  // this->initialize(vv);
}

/**
 * @brief Construct a new Phase Field Operator:: Phase Field Operator object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param spatial
 * @param params
 * @param vars
 * @param source_term_name
 */
template <class T, int DIM, class NLFI>
PhaseFieldOperatorBase<T, DIM, NLFI>::PhaseFieldOperatorBase(
    SpatialDiscretization<T, DIM> *spatial, const Parameters &params, Variables<T, DIM> &vars,
    AnalyticalFunctions<DIM> source_term_name)
    : mfem::TimeDependentOperator(spatial->getSize(), 0.0),
      M(NULL),
      N(NULL),
      reduced_oper(NULL),
      vars_(vars),
      current_dt_(0.0),
      z(height) {
  this->fespace_ = spatial->get_finite_element_space();
  this->src_func_ = source_term_name.getFunction();
}

/**
 * @brief Initialization stage (call by imeDiscretization<PST, OPE, VAR>::initialize())
 *
 * @param vv
 */
template <class T, int DIM, class NLFI>
void PhaseFieldOperatorBase<T, DIM, NLFI>::initialize(const double &initial_time) {
  Catch_Time_Section("PhaseFieldOperatorBase::initialize");
  this->SetTime(initial_time);

  auto &vv = vars_.get_variable("phi");
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
 * @brief Set current dt, unk values - needed to compute action and Jacobian.
 *solution_coef
 * @param dt time-step
 * @param u unknown vector
 */
template <class T, int DIM, class NLFI>
void PhaseFieldOperatorBase<T, DIM, NLFI>::SetTransientParameters(const double dt,
                                                                  const mfem::Vector &u) {
  Catch_Time_Section("PhaseFieldOperatorBase::SetTransientParameters");
  if (N != nullptr) {
    delete N;
  }

  if (reduced_oper != nullptr) {
    delete reduced_oper;
  }

  ////////////////////////////////////////////
  // PhaseField reduced operator N
  ////////////////////////////////////////////
  N = new mfem::NonlinearForm(this->fespace_);
  mfem::GridFunction un_gf(this->fespace_);
  un_gf.SetFromTrueDofs(u);
  mfem::GridFunction un(this->fespace_);
  un.SetFromTrueDofs(u);

  NLFI *nlfi_ptr = set_nlfi_ptr(dt, u);
  N->AddDomainIntegrator(nlfi_ptr);
  N->SetEssentialTrueDofs(this->ess_tdof_list_);

  reduced_oper = new PhaseFieldReducedOperator(M, N);

  ////////////////////////////////////////////
  // Newton Solver
  ////////////////////////////////////////////
  this->rhs_solver_ = new Solver(SolverType::NEWTON, PreconditionerType::UMFPACK, *reduced_oper);

  this->newton_solver_ = this->rhs_solver_->get_solver();
}

/**
 * @brief Set current dt, unk values - needed to compute action and Jacobian.
 *
 * @param dt time-step
 * @param u unknown vector
 * @param ess_tdof_list array of dofs
 */
template <class T, int DIM, class NLFI>
void PhaseFieldOperatorBase<T, DIM, NLFI>::SetConstantParameters(const double dt, mfem::Vector &u) {
  Catch_Time_Section("PhaseFieldOperatorBase::SetConstantParameters");
  if (M != nullptr) {
    delete M;
  }
  ////////////////
  // Mass matrix (constant)
  ////////////////
  M = new mfem::BilinearForm(this->fespace_);
  M->AddDomainIntegrator(new mfem::MassIntegrator());
  M->Assemble(0);
  M->FormSystemMatrix(this->ess_tdof_list_, Mmat);

  this->mass_solver_ = new Solver(SolverType::CG, PreconditionerType::SMOOTHER, Mmat);
  this->M_solver_ = this->mass_solver_->get_solver();
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
template <class T, int DIM, class NLFI>
void PhaseFieldOperatorBase<T, DIM, NLFI>::Mult(const mfem::Vector &u, mfem::Vector &du_dt) const {
  Catch_Time_Section("PhaseFieldOperatorBase::Mult");

  const auto sc = height;
  mfem::Vector v(u.GetData(), sc);
  mfem::Vector dv_dt(du_dt.GetData(), sc);

  N->Mult(v, z);

  // Source term
  if (!(this->src_func_ == nullptr)) {
    mfem::Vector source_term;
    this->get_source_term(source_term);
    z -= source_term;
  }

  z.Neg();  // z = -z
  this->M_solver_->Mult(z, dv_dt);
}

/**
 * @brief  Solve the Backward-Euler equation: k = f(phi + dt*k, t), for the unknown k.
 *
 * @param dt current time step
 * @param u unknown vector
 * @param du_dt unkwon time derivative vector
 */
template <class T, int DIM, class NLFI>
void PhaseFieldOperatorBase<T, DIM, NLFI>::ImplicitSolve(const double dt, const mfem::Vector &u,
                                                         mfem::Vector &du_dt) {
  Catch_Time_Section("PhaseFieldOperatorBase::ImplicitSolve");

  const auto sc = height;
  mfem::Vector v(u.GetData(), sc);
  mfem::Vector dv_dt(du_dt.GetData(), sc);
  // // Solve the equation:
  // //    du_dt = M^{-1}*[-K(u + dt*du_dt)]
  // // for du_dt

  this->bcs_->SetBoundaryConditions(v);
  this->SetTransientParameters(dt, v);

  reduced_oper->SetParameters(dt, &v);

  // Source term
  mfem::Vector source_term;
  if (!(this->src_func_ == nullptr)) {
    this->get_source_term(source_term);
  }

  dv_dt = v;
  dv_dt *= (1. / dt);
  // UtilsForDebug::memory_checkpoint("PhaseFieldOperatorBase::ImplicitSolve : before Newton Mult");
  this->newton_solver_->Mult(source_term, dv_dt);
  delete this->rhs_solver_;

  // UtilsForDebug::memory_checkpoint("PhaseFieldOperatorBase::ImplicitSolve : after Newton Mult");

  dv_dt.SetSubVector(this->ess_tdof_list_, 0.0);  // pour  Dirichlet ... uniquement?
  // std::cout << " PhaseFieldOperatorBase this->newton_solver_->Mult " << std::endl;

  MFEM_VERIFY(this->newton_solver_->GetConverged(), "Nonlinear solver did not converge.");
}

/**
 * @brief Compute L2 error
 *
 * @param u unknown vector
 * @return const double
 */
template <class T, int DIM, class NLFI>
void PhaseFieldOperatorBase<T, DIM, NLFI>::ComputeError(
    const int &it, const double &dt, const double &t, const mfem::Vector &u,
    std::function<double(const mfem::Vector &, double)> solution_func) {
  Catch_Time_Section("PhaseFieldOperatorBase::ComputeError");

  mfem::GridFunction gf(this->fespace_);
  gf.SetFromTrueDofs(u);
  mfem::FunctionCoefficient solution_coef(solution_func);
  solution_coef.SetTime(this->GetTime());
  const auto error = gf.ComputeL2Error(solution_coef);

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
void PhaseFieldOperatorBase<T, DIM, NLFI>::get_source_term(mfem::Vector &source_term) const {
  std::shared_ptr<mfem::LinearForm> RHSS(new mfem::LinearForm(this->fespace_));
  mfem::FunctionCoefficient src(this->src_func_);
  src.SetTime(this->GetTime());
  RHSS->AddDomainIntegrator(new mfem::DomainLFIntegrator(src));
  RHSS->Assemble();
  source_term = *RHSS.get();
  // source_term.Print();
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
PhaseFieldOperatorBase<T, DIM, NLFI>::get_time_specialized() const {
  return this->time_specialized_;
}

/**
 * @brief Destroy the Phase Field Operator:: Phase Field Operator object
 *
 */
template <class T, int DIM, class NLFI>
PhaseFieldOperatorBase<T, DIM, NLFI>::~PhaseFieldOperatorBase() {}

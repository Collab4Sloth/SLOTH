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
#include "Coefficients/DensityCoefficient.hpp"
#include "Coefficients/EnergyCoefficient.hpp"
#include "Coefficients/MobilityCoefficient.hpp"
#include "Coefficients/PhaseChangeFunction.hpp"
#include "Coefficients/SourceTermCoefficient.hpp"
#include "Operators/OperatorBase.hpp"
#include "Operators/ReducedOperator.hpp"
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
 * @brief PhaseFieldOperatorBase class
 *
 */
template <class T, int DIM, class NLFI>
class PhaseFieldOperatorBase : public OperatorBase<T, DIM, NLFI>,
                               public mfem::TimeDependentOperator {
 private:
  mfem::ODESolver *ode_solver_;
  bool isExplicit_{false};
  void set_ODE_solver(const TimeScheme::value &ode_solver);

  VSolverType mass_solver_;
  VSolverType mass_precond_;
  Parameters mass_solver_params_;
  Parameters mass_precond_params_;

  void set_default_mass_solver();

 protected:
  /// Mass operator

  mfem::ParGridFunction *mass_gf_;
  mfem::ParBilinearForm *M;              // mass operator
  mfem::ProductCoefficient *MassCoeff_;  // mass coefficient
  std::shared_ptr<LSolver> mass_matrix_solver_;
  VSharedMFEMSolver M_solver_;  // Krylov solver for inverting the mass matrix M

  // CCI
  //  mfem::SparseMatrix Mmat;

  mfem::HypreParMatrix *Mmat;
  // CCI
  void build_mass_matrix(const mfem::Vector &u);
  bool constant_mass_matrix_{true};

  /** Nonlinear operator defining the reduced backward Euler equation for the
      velocity. Used in the implementation of method ImplicitSolve. */
  PhaseFieldReducedOperator *reduced_oper;

 public:
  PhaseFieldOperatorBase(SpatialDiscretization<T, DIM> const *spatial, TimeScheme::value ode);

  PhaseFieldOperatorBase(SpatialDiscretization<T, DIM> const *spatial, TimeScheme::value ode,
                         AnalyticalFunctions<DIM> source_term_name);

  PhaseFieldOperatorBase(SpatialDiscretization<T, DIM> const *spatial, const Parameters &params,
                         TimeScheme::value ode);

  PhaseFieldOperatorBase(SpatialDiscretization<T, DIM> const *spatial, const Parameters &params,
                         TimeScheme::value ode, AnalyticalFunctions<DIM> source_term_name);

  void Mult(const mfem::Vector &u, mfem::Vector &du_dt) const override;
  /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
      This is the only requirement for high-order SDIRK implicit integration.*/
  void ImplicitSolve(const double dt, const mfem::Vector &u, mfem::Vector &k) override;

  virtual ~PhaseFieldOperatorBase();

  // User-defined Solvers
  void overload_mass_solver(VSolverType SOLVER, const Parameters &s_params);
  void overload_mass_solver(VSolverType SOLVER, const Parameters &s_params, VSolverType PRECOND,
                            const Parameters &p_params);
  void overload_mass_preconditioner(VSolverType PRECOND, const Parameters &p_params);

  // Virtual methods
  void set_default_properties() override = 0;

  virtual void get_mass_coefficient(const mfem::Vector &u);

  // void initialize(const double &initial_time, Variables<T, DIM> &vars) override;
  void initialize(const double &initial_time, Variables<T, DIM> &vars,
                  std::vector<Variables<T, DIM> *> auxvars) override;
  // Pure virtual methods
  void SetConstantParameters(const double dt, mfem::Vector &u) override;
  void SetTransientParameters(const double dt, const mfem::Vector &u) override;
  void solve(mfem::Vector &unk, double &next_time, const double &current_time,
             double current_time_step, const int iter) override;
  NLFI *set_nlfi_ptr(const double dt, const mfem::Vector &u) override = 0;
  void get_parameters() override = 0;
  void ComputeEnergies(const int &it, const double &t, const double &dt,
                       const mfem::Vector &u) override = 0;
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Construct a new PhaseFieldOperatorBase object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param spatial
 * @param ode
 */
template <class T, int DIM, class NLFI>
PhaseFieldOperatorBase<T, DIM, NLFI>::PhaseFieldOperatorBase(
    SpatialDiscretization<T, DIM> const *spatial, TimeScheme::value ode)
    : OperatorBase<T, DIM, NLFI>(spatial),
      mfem::TimeDependentOperator(spatial->getSize(), 0.0),
      mass_gf_(nullptr),
      M(NULL),
      MassCoeff_(NULL),
      reduced_oper(NULL) {
  this->set_ODE_solver(ode);
  this->set_default_mass_solver();
}

/**
 * @brief Construct a new PhaseFieldOperatorBase object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param spatial
 * @param ode
 * @param source_term_name
 */
template <class T, int DIM, class NLFI>
PhaseFieldOperatorBase<T, DIM, NLFI>::PhaseFieldOperatorBase(
    SpatialDiscretization<T, DIM> const *spatial, TimeScheme::value ode,
    AnalyticalFunctions<DIM> source_term_name)
    : OperatorBase<T, DIM, NLFI>(spatial, source_term_name),
      mfem::TimeDependentOperator(spatial->getSize(), 0.0),
      mass_gf_(nullptr),
      M(NULL),
      MassCoeff_(NULL),
      reduced_oper(NULL) {
  this->set_ODE_solver(ode);
  this->set_default_mass_solver();
}

/**
 * @brief Construct a new PhaseFieldOperatorBase object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param spatial
 * @param params
 * @param ode
 */
template <class T, int DIM, class NLFI>
PhaseFieldOperatorBase<T, DIM, NLFI>::PhaseFieldOperatorBase(
    SpatialDiscretization<T, DIM> const *spatial, const Parameters &params, TimeScheme::value ode)
    : OperatorBase<T, DIM, NLFI>(spatial, params),
      mfem::TimeDependentOperator(spatial->getSize(), 0.0),
      mass_gf_(nullptr),
      M(NULL),
      MassCoeff_(NULL),
      reduced_oper(NULL) {
  this->set_ODE_solver(ode);
  this->set_default_mass_solver();
}

/**
 * @brief Construct a new PhaseFieldOperatorBase object
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param spatial
 * @param params
 * @param vars
 * @param ode
 * @param auxvars
 * @param source_term_name
 */
template <class T, int DIM, class NLFI>
PhaseFieldOperatorBase<T, DIM, NLFI>::PhaseFieldOperatorBase(
    SpatialDiscretization<T, DIM> const *spatial, const Parameters &params, TimeScheme::value ode,
    AnalyticalFunctions<DIM> source_term_name)
    : OperatorBase<T, DIM, NLFI>(spatial, params, source_term_name),
      mfem::TimeDependentOperator(spatial->getSize(), 0.0),
      mass_gf_(nullptr),
      M(NULL),
      MassCoeff_(NULL),
      reduced_oper(NULL) {
  this->set_ODE_solver(ode);
  this->set_default_mass_solver();
}

/**
 * @brief  Set the ODE time marching
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param ode_solver
 */
template <class T, int DIM, class NLFI>
void PhaseFieldOperatorBase<T, DIM, NLFI>::set_ODE_solver(const TimeScheme::value &ode_solver) {
  // TODO(cci) faire un template?
  switch (ode_solver) {
    case TimeScheme::EulerExplicit: {
      this->ode_solver_ = new mfem::ForwardEulerSolver;
      this->isExplicit_ = true;
      break;
    }
    case TimeScheme::EulerImplicit: {
      this->ode_solver_ = new mfem::BackwardEulerSolver;
      this->isExplicit_ = false;
      break;
    }
    case TimeScheme::RungeKutta4: {
      this->ode_solver_ = new mfem::SDIRK33Solver;
      this->isExplicit_ = false;
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
 * @tparam NLFI
 * @param initial_time
 * @param vars
 */
template <class T, int DIM, class NLFI>
void PhaseFieldOperatorBase<T, DIM, NLFI>::initialize(const double &initial_time,
                                                      Variables<T, DIM> &vars,
                                                      std::vector<Variables<T, DIM> *> auxvars) {
  Catch_Time_Section("PhaseFieldOperatorBase::initialize");

  this->SetTime(initial_time);

  OperatorBase<T, DIM, NLFI>::initialize(initial_time, vars, auxvars);

  this->ode_solver_->Init(*this);
}

/**
 * @brief Solve the time-dependent problem
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param unk
 * @param current_time
 * @param current_time_step
 */
template <class T, int DIM, class NLFI>
void PhaseFieldOperatorBase<T, DIM, NLFI>::solve(mfem::Vector &unk, double &next_time,
                                                 const double &current_time,
                                                 double current_time_step, const int iter) {
  this->current_time_ = current_time;

  this->ode_solver_->Step(unk, next_time, current_time_step);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template <class T, int DIM, class NLFI>
void PhaseFieldOperatorBase<T, DIM, NLFI>::get_mass_coefficient(const mfem::Vector &u) {
  if (this->MassCoeff_ != nullptr) {
    delete this->MassCoeff_;
  }
  auto mass_coefficient = new mfem::ConstantCoefficient(1.0);
  this->MassCoeff_ = new mfem::ProductCoefficient(*mass_coefficient, *mass_coefficient);
}

/**
 * @brief build the mass matrix
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 */
template <class T, int DIM, class NLFI>
void PhaseFieldOperatorBase<T, DIM, NLFI>::build_mass_matrix(const mfem::Vector &u) {
  this->get_mass_coefficient(u);

  if (M != nullptr) {
    delete M;
  }
  ////////////////
  // Mass matrix (constant)
  ////////////////
  M = new mfem::ParBilinearForm(this->fespace_);

  if (!this->isExplicit_) {
    M->AddDomainIntegrator(new mfem::MassIntegrator(*this->MassCoeff_));
  } else {
    M->AddDomainIntegrator(new mfem::LumpedIntegrator(new mfem::MassIntegrator(*this->MassCoeff_)));
  }
  M->Assemble(0);
  M->Finalize(0);

  Mmat = M->ParallelAssemble();
  std::unique_ptr<mfem::HypreParMatrix> Me(Mmat->EliminateRowsCols(this->ess_tdof_list_));

  this->mass_matrix_solver_ =
      std::make_shared<LSolver>(this->mass_solver_, this->mass_solver_params_, this->mass_precond_,
                                this->mass_precond_params_, *Mmat);
  this->M_solver_ = this->mass_matrix_solver_->get_solver();
}

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

  ////////////////////////////////////////////
  // Variable mass matrix
  ////////////////////////////////////////////
  if (!this->constant_mass_matrix_) {
    this->build_mass_matrix(u);
  }
  ////////////////////////////////////////////
  // PhaseField non linear form
  ////////////////////////////////////////////
  this->build_nonlinear_form(dt, u);
  ////////////////////////////////////////////
  // PhaseField reduced operator N
  ////////////////////////////////////////////
  if (reduced_oper != nullptr) {
    delete reduced_oper;
  }
  reduced_oper = new PhaseFieldReducedOperator(M, this->N, this->ess_tdof_list_);
  ////////////////////////////////////////////
  // Newton Solver
  ////////////////////////////////////////////
  this->SetNewtonAlgorithm(reduced_oper);
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
  if (this->constant_mass_matrix_) {
    this->build_mass_matrix(u);
  }
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

  const auto sc = this->height_;
  mfem::Vector v(u.GetData(), sc);
  mfem::Vector dv_dt(du_dt.GetData(), sc);

  this->N->Mult(v, this->z);

  // Source term
  mfem::Vector source_term;
  mfem::ParLinearForm *RHS = new mfem::ParLinearForm(this->fespace_);
  if (!(this->src_func_ == nullptr)) {
    this->get_source_term(source_term, RHS);
    this->z -= source_term;
  }

  this->z.Neg();  // z = -z

  std::visit(
      [&](auto &&arg) {
        using TT = std::decay_t<decltype(arg)>;
        if constexpr (!std::is_same_v<TT, std::shared_ptr<std::monostate>>) {
          arg->Mult(this->z, dv_dt);
        }
      },
      this->M_solver_);
  delete RHS;
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

  const auto sc = this->height_;
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
  mfem::ParLinearForm *RHS = new mfem::ParLinearForm(this->fespace_);
  if (!(this->src_func_ == nullptr)) {
    this->get_source_term(source_term, RHS);
  }

  dv_dt = v;
  dv_dt *= (1. / dt);
  // UtilsForDebug::memory_checkpoint("PhaseFieldOperatorBase::ImplicitSolve : before Newton Mult");
  this->newton_solver_->Mult(source_term, dv_dt);
  delete this->rhs_solver_;
  delete RHS;

  // UtilsForDebug::memory_checkpoint("PhaseFieldOperatorBase::ImplicitSolve : after Newton
  // Mult");

  MFEM_VERIFY(this->newton_solver_->GetConverged(), "Nonlinear solver did not converge.");
}

/**
 * @brief  Overload the default options for solver used to invert the mass matrix
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param SOLVER
 * @param s_params
 */
template <class T, int DIM, class NLFI>
void PhaseFieldOperatorBase<T, DIM, NLFI>::overload_mass_solver(VSolverType SOLVER,
                                                                const Parameters &s_params) {
  this->mass_solver_ = SOLVER;
  this->mass_solver_params_ = s_params;
}

/**
 * @brief  Overload the default options for solver used to invert the mass matrix with
 * preconditionners
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param PRECOND
 * @param p_params
 */
template <class T, int DIM, class NLFI>
void PhaseFieldOperatorBase<T, DIM, NLFI>::overload_mass_solver(VSolverType SOLVER,
                                                                const Parameters &s_params,
                                                                VSolverType PRECOND,
                                                                const Parameters &p_params) {
  this->overload_solver(SOLVER, s_params);
  this->mass_precond_ = PRECOND;
  this->mass_precond_params_ = p_params;
}

/**
 * @brief  Overload the default options for preconditionners associated with the solver used to
 * invert the mass matrix
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param PRECOND
 * @param p_params
 */
template <class T, int DIM, class NLFI>
void PhaseFieldOperatorBase<T, DIM, NLFI>::overload_mass_preconditioner(
    VSolverType PRECOND, const Parameters &p_params) {
  this->mass_precond_ = PRECOND;
  this->mass_precond_params_ = p_params;
}

/**
 * @brief Set the default solver and preconditioner associated with the mass matrix
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 */
template <class T, int DIM, class NLFI>
void PhaseFieldOperatorBase<T, DIM, NLFI>::set_default_mass_solver() {
  auto s_params = Parameters(Parameter("description", "CG Solver"));
  auto p_params = Parameters(Parameter("description", "HYPRE_SMOOTHER preconditioner"));

  this->mass_solver_ = IterativeSolverType::CG;
  this->mass_solver_params_ = s_params;
  this->mass_precond_ = HyprePreconditionerType::HYPRE_SMOOTHER;
  this->mass_precond_params_ = p_params;
}
/**
 * @brief Destroy the Phase Field Operator:: Phase Field Operator object
 *
 */
template <class T, int DIM, class NLFI>
PhaseFieldOperatorBase<T, DIM, NLFI>::~PhaseFieldOperatorBase() {}

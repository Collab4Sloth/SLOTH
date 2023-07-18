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
#include <vector>
#include "BCs/BoundaryConditions.hpp"
#include "Coefficients/EnergyCoefficient.hpp"
#include "Coefficients/MobilityCoefficient.hpp"
#include "Coefficients/PhaseChangeCoefficient.hpp"
#include "Coefficients/SourceTermCoefficient.hpp"
#include "Integrators/DiffusionNLFIntegrator.hpp"
#include "Operators/ReducedOperator.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Spatial/Spatial.hpp"
#include "Utils/AnalyticalFunctions.hpp"
#include "Utils/PhaseFieldConstants.hpp"
#include "Utils/PhaseFieldOptions.hpp"
#include "Utils/UtilsForSolvers.hpp"
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
  UtilsForSolvers ut_solver_;

  // Results
  std::map<std::tuple<int, double, double>, double> error_l2_;
  std::map<std::tuple<int, double, double>, double> energy_density_;
  std::map<std::tuple<int, double, double>, double> energy_interface_;

 protected:
  mfem::FiniteElementSpace *fespace_;
  mfem::Array<int> ess_tdof_list_;

  // Mass operator
  mfem::BilinearForm *M;    // mass operator
  mfem::CGSolver M_solver;  // Krylov solver for inverting the mass matrix M
  mfem::DSmoother M_prec;   // Preconditioner for the mass matrix M
  mfem::SparseMatrix Mmat;

  // Right-Hand-Side
  mfem::LinearForm *RHS;
  mfem::NonlinearForm *N;
  mfem::Solver *J_solver;  // Solver for the Jacobian solve in the Newton method
  mfem::Solver *J_prec;    // Preconditioner for the Jacobian solve in the Newton method

  BoundaryConditions<T, DIM> *bcs_;
  Variables<T, DIM> vars_;

  mfem::NewtonSolver newton_solver_;
  /** Nonlinear operator defining the reduced backward Euler equation for the
      velocity. Used in the implementation of method ImplicitSolve. */
  PhaseFieldReducedOperator *reduced_oper;
  // mfem::Vector u_ini;
  // TODO(ci230846) : une liste de paramètres à généraliser avec la classe Parameters
  double mobility_coeff_, omega_, lambda_;

  std::function<double(const mfem::Vector &, double)> src_func_;

  double current_dt_;
  mutable mfem::Vector z;  // auxiliary vector
  std::string source_term_name_;

 public:
  PhaseFieldOperatorBase(SpatialDiscretization<T, DIM> *spatial, const Parameters &params,
                         Variables<T, DIM> &vars);

  template <typename... Args>
  PhaseFieldOperatorBase(SpatialDiscretization<T, DIM> *spatial, const Parameters &params,
                         Variables<T, DIM> &vars, const std::string &source_term_name,
                         Args... args);

  virtual void Mult(const mfem::Vector &u, mfem::Vector &du_dt) const;
  /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
      This is the only requirement for high-order SDIRK implicit integration.*/
  virtual void ImplicitSolve(const double dt, const mfem::Vector &u, mfem::Vector &k);
  void SetConstantParameters(const double dt, mfem::Vector &u);
  void SetTransientParameters(const double dt, const mfem::Vector &u);
  void initialize(const double &initial_time);
  void ComputeEnergies(const std::tuple<int, double, double> &iter, const mfem::Vector &u);
  void ComputeError(const std::tuple<int, double, double> &iter, const mfem::Vector &u,
                    std::function<double(const mfem::Vector &, double)> solution_func);
  void get_source_term(mfem::Vector &source_term) const;

  const std::map<std::tuple<int, double, double>, double> get_l2_error() const;
  const std::map<std::tuple<int, double, double>, double> get_energy_density() const;
  const std::map<std::tuple<int, double, double>, double> get_energy_interface() const;

  virtual ~PhaseFieldOperatorBase();

  // Pure virtual methods
  virtual NLFI *set_nlfi_ptr(const double dt, const mfem::Vector &u) = 0;
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/**
 * @brief Construct a new Phase Field Operator:: Phase Field Operator object
 *
 * @param fespace Finite Element space
 * @param params list of Parameters
 * @param u unknown vector
 */
template <class T, int DIM, class NLFI>
PhaseFieldOperatorBase<T, DIM, NLFI>::PhaseFieldOperatorBase(SpatialDiscretization<T, DIM> *spatial,
                                                             const Parameters &params,
                                                             Variables<T, DIM> &vars)
    : mfem::TimeDependentOperator(spatial->getSize(), 0.0),
      M(NULL),
      N(NULL),
      vars_(vars),
      current_dt_(0.0),
      z(height) {
  this->fespace_ = spatial->get_finite_element_space();
  this->omega_ = params.get_parameter_value("omega");
  this->lambda_ = params.get_parameter_value("lambda");
  this->mobility_coeff_ = params.get_parameter_value("mobility");
}
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/**
 * @brief Construct a new Phase Field Operator:: Phase Field Operator object
 *
 * @param fespace Finite Element space
 * @param params list of Parameters
 * @param u unknown vector
 */
template <class T, int DIM, class NLFI>
template <typename... Args>
PhaseFieldOperatorBase<T, DIM, NLFI>::PhaseFieldOperatorBase(SpatialDiscretization<T, DIM> *spatial,
                                                             const Parameters &params,
                                                             Variables<T, DIM> &vars,
                                                             const std::string &source_term_name,
                                                             Args... args)
    : mfem::TimeDependentOperator(spatial->getSize(), 0.0),
      M(NULL),
      N(NULL),
      vars_(vars),
      current_dt_(0.0),
      z(height),
      source_term_name_(source_term_name) {
  this->fespace_ = spatial->get_finite_element_space();
  this->omega_ = params.get_parameter_value("omega");
  this->lambda_ = params.get_parameter_value("lambda");
  this->mobility_coeff_ = params.get_parameter_value("mobility");
  AnalyticalFunctions<DIM> source_func;
  this->src_func_ = source_func.getAnalyticalFunctions(this->source_term_name_, args...);

  // auto &vv = vars.get_variable("phi");
  // this->initialize(vv);
}

/**
 * @brief Initialization stage (call by imeDiscretization<PST, OPE, VAR>::initialize())
 *
 * @param vv
 */
template <class T, int DIM, class NLFI>
void PhaseFieldOperatorBase<T, DIM, NLFI>::initialize(const double &initial_time) {
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
  delete N;
  delete reduced_oper;

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
  this->ut_solver_.SetSolverParameters(
      this->newton_solver_, NewtonDefaultConstant::print_level,
      NewtonDefaultConstant::iterative_mode, NewtonDefaultConstant::iter_max,
      NewtonDefaultConstant::rel_tol, NewtonDefaultConstant::abs_tol);
  // TODO(ci230846) : cette partie devra etre generalisee pour un solveur iteratif
  J_solver = new mfem::UMFPackSolver;
  this->ut_solver_.BuildSolver(this->newton_solver_, *J_solver, *reduced_oper);
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
  delete M;
  ////////////////
  // Mass matrix (constant)
  ////////////////
  M = new mfem::BilinearForm(this->fespace_);
  M->AddDomainIntegrator(new mfem::MassIntegrator());
  M->Assemble(0);
  mfem::SparseMatrix tmp;
  M->FormSystemMatrix(this->ess_tdof_list_, Mmat);
  this->ut_solver_.SetSolverParameters(
      M_solver, MassDefaultConstant::print_level, MassDefaultConstant::iterative_mode,
      MassDefaultConstant::iter_max, MassDefaultConstant::rel_tol, MassDefaultConstant::abs_tol);
  this->ut_solver_.BuildSolver(M_solver, M_prec, Mmat);
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
  const auto sc = height;
  mfem::Vector v(u.GetData(), sc);
  mfem::Vector dv_dt(du_dt.GetData(), sc);

  N->Mult(v, z);

  // Source term
  if (!this->source_term_name_.empty()) {
    mfem::Vector source_term;
    this->get_source_term(source_term);
    z -= source_term;
  }

  z.Neg();  // z = -z
  M_solver.Mult(z, dv_dt);
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
  if (!this->source_term_name_.empty()) {
    this->get_source_term(source_term);
  }

  dv_dt = v;
  dv_dt *= (1. / dt);
  this->newton_solver_.Mult(source_term, dv_dt);
  dv_dt.SetSubVector(this->ess_tdof_list_, 0.0);  // pour  Dirichlet ... uniquement?
  // std::cout << " PhaseFieldOperatorBase this->newton_solver_->Mult " << std::endl;

  MFEM_VERIFY(this->newton_solver_.GetConverged(), "Nonlinear solver did not converge.");
}

/**
 * @brief Compute Phase-field Energies
 *
 * @param u unknown vector
 * @return const double
 */
template <class T, int DIM, class NLFI>
void PhaseFieldOperatorBase<T, DIM, NLFI>::ComputeEnergies(
    const std::tuple<int, double, double> &iter, const mfem::Vector &u) {
  mfem::GridFunction un_gf(this->fespace_);
  un_gf.SetFromTrueDofs(u);
  mfem::GridFunction gf(this->fespace_);
  EnergyCoefficient g(&un_gf, 0.5 * this->lambda_, this->omega_);
  gf.ProjectCoefficient(g);
  mfem::GridFunction sigf(this->fespace_);
  EnergyCoefficient sig(&un_gf, this->lambda_, 0.);
  sigf.ProjectCoefficient(sig);

  // Calcul de l'intégrale de l'objet FunctionCoefficient sur le domaine
  mfem::ConstantCoefficient zero(0.);
  const auto energy = gf.ComputeL1Error(zero);
  this->energy_density_.try_emplace(iter, energy);
  const auto interfacial_energy = sigf.ComputeL1Error(zero);
  this->energy_interface_.try_emplace(iter, interfacial_energy);
}

/**
 * @brief Compute L2 error
 *
 * @param u unknown vector
 * @return const double
 */
template <class T, int DIM, class NLFI>
void PhaseFieldOperatorBase<T, DIM, NLFI>::ComputeError(
    const std::tuple<int, double, double> &iter, const mfem::Vector &u,
    std::function<double(const mfem::Vector &, double)> solution_func) {
  mfem::GridFunction gf(this->fespace_);
  gf.SetFromTrueDofs(u);
  mfem::FunctionCoefficient solution_coef(solution_func);
  solution_coef.SetTime(this->GetTime());
  const auto error = gf.ComputeL2Error(solution_coef);
  this->error_l2_.try_emplace(iter, error);
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
}

/**
 * @brief Get a value of a map of L2 errors at a given iteration
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param iter
 * @return const std::map<std::tuple<int, double, double>, double>
 */
template <class T, int DIM, class NLFI>
const std::map<std::tuple<int, double, double>, double>
PhaseFieldOperatorBase<T, DIM, NLFI>::get_l2_error() const {
  return this->error_l2_;
}

/**
 * @brief Get a of a map of interfacial energy  at a given iteration
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @return const std::map<std::tuple<int, double, double>, double>
 */
template <class T, int DIM, class NLFI>
const std::map<std::tuple<int, double, double>, double>
PhaseFieldOperatorBase<T, DIM, NLFI>::get_energy_density() const {
  return this->energy_density_;
}

/**
 * @brief Get a of a map of interfacial energy  at a given iteration
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @return const std::map<std::tuple<int, double, double>, double>
 */
template <class T, int DIM, class NLFI>
const std::map<std::tuple<int, double, double>, double>
PhaseFieldOperatorBase<T, DIM, NLFI>::get_energy_interface() const {
  return this->energy_interface_;
}

/**
 * @brief Destroy the Phase Field Operator:: Phase Field Operator object
 *
 */
template <class T, int DIM, class NLFI>
PhaseFieldOperatorBase<T, DIM, NLFI>::~PhaseFieldOperatorBase() {
  // delete J_solver;
  // delete J_prec;
  // delete reduced_oper;
}

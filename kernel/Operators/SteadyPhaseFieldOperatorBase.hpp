/**
 * @file SteadyPhaseFieldOperatorBase.hpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief Steady PhaseField Operator base
 * @version 0.1
 * @date 2024-07-19
 *
 * @copyright Copyright (c) 2024
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
#include "Operators/OperatorBase.hpp"
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
 * @brief SteadyPhaseFieldOperatorBase class
 *
 */
template <class T, int DIM, class NLFI>
class SteadyPhaseFieldOperatorBase : public OperatorBase<T, DIM, NLFI> {
 public:
  explicit SteadyPhaseFieldOperatorBase(SpatialDiscretization<T, DIM> const *spatial);

  SteadyPhaseFieldOperatorBase(SpatialDiscretization<T, DIM> const *spatial,
                               AnalyticalFunctions<DIM> source_term_name);

  SteadyPhaseFieldOperatorBase(SpatialDiscretization<T, DIM> const *spatial,
                               const Parameters &params);

  SteadyPhaseFieldOperatorBase(SpatialDiscretization<T, DIM> const *spatial,
                               const Parameters &params, AnalyticalFunctions<DIM> source_term_name);
  /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
      This is the only requirement for high-order SDIRK implicit integration.*/

  void Mult(const mfem::Vector &k, mfem::Vector &y) const override;
  mfem::Operator &GetGradient(const mfem::Vector &k) const override;

  virtual ~SteadyPhaseFieldOperatorBase();

  // Virtual methods
  void set_default_properties() override = 0;

  // void initialize(const double &initial_time, Variables<T, DIM> &vars) override;
  void initialize(const double &initial_time, Variables<T, DIM> &vars,
                  std::vector<Variables<T, DIM> *> auxvars) override;
  // Pure virtual methods
  void SetConstantParameters(const double dt, mfem::Vector &u) override;
  void SetTransientParameters(const double dt, const mfem::Vector &u) override;
  void solve(mfem::Vector &unk, double &next_time, const double &current_time, double dt,
             const int iter) override;
  NLFI *set_nlfi_ptr(const double dt, const mfem::Vector &u) override = 0;
  void get_parameters() override = 0;
  void ComputeEnergies(const int &it, const double &t, const double &dt,
                       const mfem::Vector &u) override = 0;
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Construct a new SteadyPhaseFieldOperatorBase object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param spatial
 */
template <class T, int DIM, class NLFI>
SteadyPhaseFieldOperatorBase<T, DIM, NLFI>::SteadyPhaseFieldOperatorBase(
    SpatialDiscretization<T, DIM> const *spatial)
    : OperatorBase<T, DIM, NLFI>(spatial) {
  const Parameters nl_parameters =
      Parameters(Parameter("description", "Newton Algorithm"), Parameter("iterative_mode", true));
  this->overload_nl_solver(NLSolverType::NEWTON, nl_parameters);
}

/**
 * @brief Construct a new SteadyPhaseFieldOperatorBase object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param spatial
 * @param source_term_name
 */
template <class T, int DIM, class NLFI>
SteadyPhaseFieldOperatorBase<T, DIM, NLFI>::SteadyPhaseFieldOperatorBase(
    SpatialDiscretization<T, DIM> const *spatial, AnalyticalFunctions<DIM> source_term_name)
    : OperatorBase<T, DIM, NLFI>(spatial, source_term_name) {
  const Parameters nl_parameters =
      Parameters(Parameter("description", "Newton Algorithm"), Parameter("iterative_mode", true));
  this->overload_nl_solver(NLSolverType::NEWTON, nl_parameters);
}

/**
 * @brief Construct a new SteadyPhaseFieldOperatorBase object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param spatial
 * @param params
 */
template <class T, int DIM, class NLFI>
SteadyPhaseFieldOperatorBase<T, DIM, NLFI>::SteadyPhaseFieldOperatorBase(
    SpatialDiscretization<T, DIM> const *spatial, const Parameters &params)
    : OperatorBase<T, DIM, NLFI>(spatial, params) {
  const Parameters nl_parameters =
      Parameters(Parameter("description", "Newton Algorithm"), Parameter("iterative_mode", true));
  this->overload_nl_solver(NLSolverType::NEWTON, nl_parameters);
}

/**
 * @brief Construct a new SteadyPhaseFieldOperatorBase object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param spatial
 * @param params
 * @param source_term_name
 */
template <class T, int DIM, class NLFI>
SteadyPhaseFieldOperatorBase<T, DIM, NLFI>::SteadyPhaseFieldOperatorBase(
    SpatialDiscretization<T, DIM> const *spatial, const Parameters &params,
    AnalyticalFunctions<DIM> source_term_name)
    : OperatorBase<T, DIM, NLFI>(spatial, params, source_term_name) {
  const Parameters nl_parameters =
      Parameters(Parameter("description", "Newton Algorithm"), Parameter("iterative_mode", true));
  this->overload_nl_solver(NLSolverType::NEWTON, nl_parameters);
}

/**
 * @brief Initialization stage (call by imeDiscretization<PST, OPE, VAR>::initialize())
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param initial_time
 * @param vars
 */
template <class T, int DIM, class NLFI>
void SteadyPhaseFieldOperatorBase<T, DIM, NLFI>::initialize(
    const double &initial_time, Variables<T, DIM> &vars, std::vector<Variables<T, DIM> *> auxvars) {
  Catch_Time_Section("SteadyPhaseFieldOperatorBase::initialize");

  OperatorBase<T, DIM, NLFI>::initialize(initial_time, vars, auxvars);
}

/**
 * @brief Solve the steady problem
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param unk
 * @param current_time
 * @param dt
 */
template <class T, int DIM, class NLFI>
void SteadyPhaseFieldOperatorBase<T, DIM, NLFI>::solve(mfem::Vector &unk, double &next_time,
                                                       const double &current_time, double dt,
                                                       const int iter) {
  this->current_time_ = current_time;
  // height (defined mfem::Operator class)
  const auto sc = this->height_;
  mfem::Vector v(unk.GetData(), sc);

  this->bcs_->SetBoundaryConditions(v);
  this->SetTransientParameters(dt, v);

  // Source term
  mfem::Vector source_term;
  mfem::ParLinearForm *RHS = new mfem::ParLinearForm(this->fespace_);
  if (!(this->src_func_ == nullptr)) {
    this->get_source_term(source_term, RHS);
  }

  this->newton_solver_->Mult(source_term, v);
  delete this->rhs_solver_;
  delete RHS;

  // current_time += dt;

  // this->current_time_ = current_time + static_cast<double>(iter) * dt;
  next_time = current_time + static_cast<double>(iter) * dt;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief Set current dt, unk values - needed to compute action and Jacobian.
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param dt
 * @param u
 */
template <class T, int DIM, class NLFI>
void SteadyPhaseFieldOperatorBase<T, DIM, NLFI>::SetConstantParameters(const double dt,
                                                                       mfem::Vector &u) {
  // Catch_Time_Section("OperatorBase::SetConstantParameters");
  // Nothing to be done
}
/**
 * @brief Set current dt, unk values - needed to compute action and Jacobian.
 *solution_coef
 * @param dt time-step
 * @param u unknown vector
 */
template <class T, int DIM, class NLFI>
void SteadyPhaseFieldOperatorBase<T, DIM, NLFI>::SetTransientParameters(const double dt,
                                                                        const mfem::Vector &u) {
  Catch_Time_Section("SteadyPhaseFieldOperatorBase::SetTransientParameters");
  ////////////////////////////////////////////
  // PhaseField non linear form
  ////////////////////////////////////////////
  this->build_nonlinear_form(dt, u);
  ////////////////////////////////////////////
  // Newton Solver
  ////////////////////////////////////////////
  this->SetNewtonAlgorithm(this->N);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief Compute RHS
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param k
 * @param y
 */
template <class T, int DIM, class NLFI>
void SteadyPhaseFieldOperatorBase<T, DIM, NLFI>::Mult(const mfem::Vector &k,
                                                      mfem::Vector &y) const {
  this->N->Mult(k, y);
  y.SetSubVector(this->ess_tdof_list_, 0.0);
}

/**
 * @brief Compute Jacobian
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param k
 * @return mfem::Operator&
 */
template <class T, int DIM, class NLFI>
mfem::Operator &SteadyPhaseFieldOperatorBase<T, DIM, NLFI>::GetGradient(
    const mfem::Vector &k) const {
  // TODO(cci) : pourquoi pas GetLocalGradient? pas de eleminaterowcols
  // mfem::HypreParMatrix *localJ =
  //     dynamic_cast<mfem::HypreParMatrix *>(&this->N->GetLocalGradient(k));

  return this->N->GetGradient(k);
}

/**
 * @brief Destroy the Steady Phase Field Operator Base< T,  DIM,  NLFI>:: Steady Phase Field
 * Operator Base object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 */
template <class T, int DIM, class NLFI>
SteadyPhaseFieldOperatorBase<T, DIM, NLFI>::~SteadyPhaseFieldOperatorBase() {}

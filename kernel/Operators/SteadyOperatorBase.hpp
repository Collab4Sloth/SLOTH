/**
 * @file SteadyOperatorBase.hpp
 * @author ci230846 (clement.introini@cea.fr)
 * @brief Steady PhaseField Operator base
 * @version 0.1
 * @date 2024-07-19
 *
 * @copyright Copyright (c) 2024
 *
 */

#include <memory>
#include <vector>

#include "AnalyticalFunctions/AnalyticalFunctions.hpp"
#include "Operators/OperatorBase.hpp"
#include "Operators/SteadyReducedOperator.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameters.hpp"
#include "Profiling/Profiling.hpp"
#include "Spatial/Spatial.hpp"
#include "Utils/Utils.hpp"
#include "Variables/Variables.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

/**
 * @brief SteadyOperatorBase class
 *
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
class SteadyOperatorBase : public OperatorBase<T, DIM, NLFI, LHS_NLFI> {
 private:
  SteadyPhaseFieldReducedOperator* steady_reduced_oper;

 public:
  explicit SteadyOperatorBase(std::vector<SpatialDiscretization<T, DIM>*> spatials);

  SteadyOperatorBase(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                     std::vector<AnalyticalFunctions<DIM>> source_term_name);

  SteadyOperatorBase(std::vector<SpatialDiscretization<T, DIM>*> spatials,
                     const Parameters& params);

  SteadyOperatorBase(std::vector<SpatialDiscretization<T, DIM>*> spatials, const Parameters& params,
                     std::vector<AnalyticalFunctions<DIM>> source_term_name);

  void Mult(const mfem::Vector& k, mfem::Vector& y) const override;
  // mfem::Operator &GetGradient(const mfem::Vector &k) const override;

  virtual ~SteadyOperatorBase();

  // Virtual methods
  void set_default_properties() override = 0;

  // void initialize(const double &initial_time, Variables<T, DIM> &vars) override;
  void initialize(const double& initial_time, Variables<T, DIM>& vars,
                  std::vector<Variables<T, DIM>*> auxvars) override;
  // Pure virtual methods
  void SetConstantParameters(const double dt, const std::vector<mfem::Vector>& u_vect) override;
  void SetTransientParameters(const double dt, const std::vector<mfem::Vector>& u_vet) override;
  void solve(std::vector<std::unique_ptr<mfem::Vector>>& vect_unk, double& next_time,
             const double& current_time, double dt, const int iter) override;
  NLFI* set_nlfi_ptr(const double dt, const std::vector<mfem::Vector>& u) override = 0;
  void get_parameters() override = 0;
  void ComputeEnergies(const int& it, const double& t, const double& dt,
                       const std::vector<mfem::Vector>& u) override = 0;
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Construct a new SteadyOperatorBase object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param spatial
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
SteadyOperatorBase<T, DIM, NLFI, LHS_NLFI>::SteadyOperatorBase(
    std::vector<SpatialDiscretization<T, DIM>*> spatials)
    : OperatorBase<T, DIM, NLFI, LHS_NLFI>(spatials), steady_reduced_oper(NULL) {
  const Parameters nl_parameters =
      Parameters(Parameter("description", "Newton Algorithm"), Parameter("iterative_mode", true));
  this->overload_nl_solver(NLSolverType::NEWTON, nl_parameters);
}

/**
 * @brief Construct a new SteadyOperatorBase object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param spatial
 * @param source_term_name
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
SteadyOperatorBase<T, DIM, NLFI, LHS_NLFI>::SteadyOperatorBase(
    std::vector<SpatialDiscretization<T, DIM>*> spatials,
    std::vector<AnalyticalFunctions<DIM>> source_term_name)
    : OperatorBase<T, DIM, NLFI, LHS_NLFI>(spatials, source_term_name), steady_reduced_oper(NULL) {
  const Parameters nl_parameters =
      Parameters(Parameter("description", "Newton Algorithm"), Parameter("iterative_mode", true));
  this->overload_nl_solver(NLSolverType::NEWTON, nl_parameters);
}

/**
 * @brief Construct a new SteadyOperatorBase object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param spatial
 * @param params
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
SteadyOperatorBase<T, DIM, NLFI, LHS_NLFI>::SteadyOperatorBase(
    std::vector<SpatialDiscretization<T, DIM>*> spatials, const Parameters& params)
    : OperatorBase<T, DIM, NLFI, LHS_NLFI>(spatials, params), steady_reduced_oper(NULL) {
  const Parameters nl_parameters =
      Parameters(Parameter("description", "Newton Algorithm"), Parameter("iterative_mode", true));
  this->overload_nl_solver(NLSolverType::NEWTON, nl_parameters);
}

/**
 * @brief Construct a new SteadyOperatorBase object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param spatial
 * @param params
 * @param source_term_name
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
SteadyOperatorBase<T, DIM, NLFI, LHS_NLFI>::SteadyOperatorBase(
    std::vector<SpatialDiscretization<T, DIM>*> spatials, const Parameters& params,
    std::vector<AnalyticalFunctions<DIM>> source_term_name)
    : OperatorBase<T, DIM, NLFI, LHS_NLFI>(spatials, params, source_term_name),
      steady_reduced_oper(NULL) {
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
template <class T, int DIM, class NLFI, class LHS_NLFI>
void SteadyOperatorBase<T, DIM, NLFI, LHS_NLFI>::initialize(
    const double& initial_time, Variables<T, DIM>& vars, std::vector<Variables<T, DIM>*> auxvars) {
  Catch_Time_Section("SteadyOperatorBase::initialize");

  OperatorBase<T, DIM, NLFI, LHS_NLFI>::initialize(initial_time, vars, auxvars);
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
template <class T, int DIM, class NLFI, class LHS_NLFI>
void SteadyOperatorBase<T, DIM, NLFI, LHS_NLFI>::solve(
    std::vector<std::unique_ptr<mfem::Vector>>& vect_unk, double& next_time,
    const double& current_time, double dt, const int iter) {
  this->current_time_ = current_time;
  this->current_dt_ = dt;

  //// Constructing array of offsets
  const size_t unk_size = vect_unk.size();

  //// Constructing BlockVector
  mfem::BlockVector block_unk(this->block_trueOffsets_);
  std::vector<mfem::Vector> u_vect;
  for (size_t i = 0; i < unk_size; i++) {
    auto& unk_i = *(vect_unk[i]);
    mfem::Vector& bb = block_unk.GetBlock(i);
    bb = unk_i;
    u_vect.emplace_back(unk_i);
  }

  this->SetTransientParameters(dt, u_vect);

  /// Apply BCs: check if need to be uncomment
  // for (size_t i = 0; i < unk_size; i++) {
  //   auto &unk_i = *(vect_unk[i]);
  //   this->bcs_[i]->SetBoundaryConditions(unk_i);
  // }

  // Source term
  // const mfem::Array<int> offsets = this->RHS->GetBlockOffsets();
  const int fes_size = this->block_trueOffsets_.Size() - 1;
  mfem::BlockVector source_term(this->block_trueOffsets_);
  source_term = 0.0;
  if (!this->src_func_.empty()) {
    for (int i = 0; i < fes_size; ++i) {
      if (this->src_func_[i] != nullptr) {
        mfem::ParLinearForm* RHS = new mfem::ParLinearForm(this->fes_[i]);
        mfem::Vector& src_i = source_term.GetBlock(i);
        this->get_source_term(i, this->src_func_[i], src_i, RHS);
        delete RHS;
      }
    }
  }
  this->newton_solver_->Mult(source_term, block_unk);
  delete this->rhs_solver_;

  for (size_t i = 0; i < unk_size; i++) {
    auto& unk_i = *(vect_unk[i]);
    const mfem::Vector& bb = block_unk.GetBlock(i);
    unk_i = bb;
  }
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
template <class T, int DIM, class NLFI, class LHS_NLFI>
void SteadyOperatorBase<T, DIM, NLFI, LHS_NLFI>::SetConstantParameters(
    const double dt, const std::vector<mfem::Vector>& u_vect) {
  // Catch_Time_Section("OperatorBase::SetConstantParameters");
  // Nothing to be done
}
/**
 * @brief Set current dt, unk values - needed to compute action and Jacobian.
 *solution_coef
 * @param dt time-step
 * @param u unknown vector
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void SteadyOperatorBase<T, DIM, NLFI, LHS_NLFI>::SetTransientParameters(
    const double dt, const std::vector<mfem::Vector>& u_vect) {
  Catch_Time_Section("SteadyOperatorBase::SetTransientParameters");

  ////////////////////////////////////////////
  //  Build the RHS of the PDEs
  ////////////////////////////////////////////
  this->build_rhs_nonlinear_form(dt, u_vect);

  ////////////////////////////////////////////
  // Build Newton Linear system
  ////////////////////////////////////////////
  if (steady_reduced_oper != nullptr) {
    delete steady_reduced_oper;
  }
  steady_reduced_oper = new SteadyPhaseFieldReducedOperator(this->RHS, this->ess_tdof_list_);

  ////////////////////////////////////////////
  // Newton Solver
  ////////////////////////////////////////////
  this->SetNewtonAlgorithm(steady_reduced_oper);
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
template <class T, int DIM, class NLFI, class LHS_NLFI>
void SteadyOperatorBase<T, DIM, NLFI, LHS_NLFI>::Mult(const mfem::Vector& k,
                                                      mfem::Vector& y) const {
  // Nothing to be done because of manage by steadyreducedoperator
}

/**
 * @brief Destroy the Steady Phase Field Operator Base< T,  DIM,  NLFI>:: Steady Phase Field
 * Operator Base object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
SteadyOperatorBase<T, DIM, NLFI, LHS_NLFI>::~SteadyOperatorBase() {}

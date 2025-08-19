/**
 * @file TransientOperatorBase.hpp
 * @author ci230846
 * @brief
 * @version 0.1
 * @date 2025-07-09
 *
 * @copyright Copyright © CEA 2025
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
#include "Coefficients/DensityCoefficient.hpp"
#include "Coefficients/EnergyCoefficient.hpp"
#include "Coefficients/MobilityCoefficient.hpp"
#include "Coefficients/PhaseChangeFunction.hpp"
#include "Operators/OperatorBase.hpp"
#include "Operators/ReducedOperator.hpp"
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
 * @brief TransientOperatorBase class
 *
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
class TransientOperatorBase : public OperatorBase<T, DIM, NLFI, LHS_NLFI>,
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
  std::vector<VSharedMFEMSolver> M_solver_;  // Krylov solver for inverting the mass matrix )M
  LHS_NLFI *lhs_nlfi_ptr_;

  /// Left-Hand-Side
  mfem::ParBlockNonlinearForm *LHS;

  // CCI
  //  mfem::SparseMatrix Mmat;
  mfem::HypreParMatrix *Mmat;
  // CCI
  void build_mass_matrix(const std::vector<mfem::Vector> &u_vect);
  void build_lhs_nonlinear_form(const double dt, const std::vector<mfem::Vector> &u);
  bool constant_mass_matrix_{true};

  // Reduced Operator
  PhaseFieldReducedOperator *reduced_oper;

 public:
  TransientOperatorBase(std::vector<SpatialDiscretization<T, DIM> *> spatials,
                        TimeScheme::value ode);

  TransientOperatorBase(std::vector<SpatialDiscretization<T, DIM> *> spatials,
                        TimeScheme::value ode,
                        std::vector<AnalyticalFunctions<DIM>> source_term_name);

  TransientOperatorBase(std::vector<SpatialDiscretization<T, DIM> *> spatials,
                        const Parameters &params, TimeScheme::value ode);

  TransientOperatorBase(std::vector<SpatialDiscretization<T, DIM> *> spatials,
                        const Parameters &params, TimeScheme::value ode,
                        std::vector<AnalyticalFunctions<DIM>> source_term_name);

  void Mult(const mfem::Vector &u, mfem::Vector &du_dt) const override;
  void ImplicitSolve(const double dt, const mfem::Vector &u, mfem::Vector &k) override;

  virtual ~TransientOperatorBase();

  // User-defined Solvers
  void overload_mass_solver(VSolverType SOLVER);
  void overload_mass_solver(VSolverType SOLVER, const Parameters &s_params);
  void overload_mass_preconditioner(VSolverType PRECOND);
  void overload_mass_preconditioner(VSolverType PRECOND, const Parameters &p_params);

  void SetExplicitTransientParameters(const std::vector<mfem::Vector> &un_vect);

  // Virtual methods
  void set_default_properties() override = 0;

  virtual void get_mass_coefficient(const mfem::Vector &u);

  // void initialize(const double &initial_time, Variables<T, DIM> &vars) override;
  void initialize(const double &initial_time, Variables<T, DIM> &vars,
                  std::vector<Variables<T, DIM> *> auxvars) override;
  // Pure virtual methods
  void SetConstantParameters(const double dt, const std::vector<mfem::Vector> &u_vect) override;
  void SetTransientParameters(const double dt, const std::vector<mfem::Vector> &u_vect) override;
  void solve(std::vector<std::unique_ptr<mfem::Vector>> &vect_unk, double &next_time,
             const double &current_time, double current_time_step, const int iter) override;
  NLFI *set_nlfi_ptr(const double dt, const std::vector<mfem::Vector> &u) override = 0;
  LHS_NLFI *set_lhs_nlfi_ptr(const double dt, const std::vector<mfem::Vector> &u);

  void get_parameters() override = 0;
  void ComputeEnergies(const int &it, const double &t, const double &dt,
                       const std::vector<mfem::Vector> &u) override = 0;
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Construct a new TransientOperatorBase object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param spatial
 * @param ode
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
TransientOperatorBase<T, DIM, NLFI, LHS_NLFI>::TransientOperatorBase(
    std::vector<SpatialDiscretization<T, DIM> *> spatials, TimeScheme::value ode)
    : OperatorBase<T, DIM, NLFI, LHS_NLFI>(spatials),
      mfem::TimeDependentOperator(this->compute_total_height(spatials),
                                  this->compute_total_width(spatials), 0.0),
      mass_gf_(nullptr),
      M(NULL),
      MassCoeff_(NULL),
      LHS(NULL),
      lhs_nlfi_ptr_(nullptr),
      reduced_oper(NULL) {
  this->set_ODE_solver(ode);
  this->set_default_mass_solver();
}

/**
 * @brief Construct a new TransientOperatorBase object
 * object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @tparam LHS_NLFI
 * @param spatials
 * @param ode
 * @param source_term_name
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
TransientOperatorBase<T, DIM, NLFI, LHS_NLFI>::TransientOperatorBase(
    std::vector<SpatialDiscretization<T, DIM> *> spatials, TimeScheme::value ode,
    std::vector<AnalyticalFunctions<DIM>> source_term_name)
    : OperatorBase<T, DIM, NLFI, LHS_NLFI>(spatials, source_term_name),
      mfem::TimeDependentOperator(this->compute_total_height(spatials),
                                  this->compute_total_width(spatials), 0.0),
      mass_gf_(nullptr),
      M(NULL),
      MassCoeff_(NULL),
      LHS(NULL),
      lhs_nlfi_ptr_(nullptr),
      reduced_oper(NULL) {
  this->set_ODE_solver(ode);
  this->set_default_mass_solver();
}

/**
 * @brief Construct a new TransientOperatorBase object
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param spatial
 * @param params
 * @param ode
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
TransientOperatorBase<T, DIM, NLFI, LHS_NLFI>::TransientOperatorBase(
    std::vector<SpatialDiscretization<T, DIM> *> spatials, const Parameters &params,
    TimeScheme::value ode)
    : OperatorBase<T, DIM, NLFI, LHS_NLFI>(spatials, params),
      mfem::TimeDependentOperator(this->compute_total_height(spatials),
                                  this->compute_total_width(spatials), 0.0),
      mass_gf_(nullptr),
      M(NULL),
      MassCoeff_(NULL),
      LHS(NULL),
      lhs_nlfi_ptr_(nullptr),
      reduced_oper(NULL) {
  this->set_ODE_solver(ode);
  this->set_default_mass_solver();
}

/**
 * @brief Construct a new TransientOperatorBase object
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
template <class T, int DIM, class NLFI, class LHS_NLFI>
TransientOperatorBase<T, DIM, NLFI, LHS_NLFI>::TransientOperatorBase(
    std::vector<SpatialDiscretization<T, DIM> *> spatials, const Parameters &params,
    TimeScheme::value ode, std::vector<AnalyticalFunctions<DIM>> source_term_name)
    : OperatorBase<T, DIM, NLFI, LHS_NLFI>(spatials, params, source_term_name),
      mfem::TimeDependentOperator(this->compute_total_height(spatials),
                                  this->compute_total_width(spatials), 0.0),
      mass_gf_(nullptr),
      M(NULL),
      LHS(NULL),
      MassCoeff_(NULL),
      reduced_oper(NULL),
      lhs_nlfi_ptr_(nullptr) {
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
template <class T, int DIM, class NLFI, class LHS_NLFI>
void TransientOperatorBase<T, DIM, NLFI, LHS_NLFI>::set_ODE_solver(
    const TimeScheme::value &ode_solver) {
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
template <class T, int DIM, class NLFI, class LHS_NLFI>
void TransientOperatorBase<T, DIM, NLFI, LHS_NLFI>::initialize(
    const double &initial_time, Variables<T, DIM> &vars, std::vector<Variables<T, DIM> *> auxvars) {
  Catch_Time_Section("TransientOperatorBase::initialize");

  this->SetTime(initial_time);

  OperatorBase<T, DIM, NLFI, LHS_NLFI>::initialize(initial_time, vars, auxvars);

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
template <class T, int DIM, class NLFI, class LHS_NLFI>
void TransientOperatorBase<T, DIM, NLFI, LHS_NLFI>::solve(
    std::vector<std::unique_ptr<mfem::Vector>> &vect_unk, double &next_time,
    const double &current_time, double current_time_step, const int iter) {
  //// Constructing array of offsets
  const size_t unk_size = vect_unk.size();

  //// Constructing BlockVector
  mfem::BlockVector block_unk(this->block_trueOffsets_);
  for (size_t i = 0; i < unk_size; i++) {
    auto &unk_i = *(vect_unk[i]);
    mfem::Vector &bb = block_unk.GetBlock(i);
    bb = unk_i;
  }
  //// Call ODE solver
  this->current_time_ = current_time;
  this->current_dt_ = current_time_step;
  this->ode_solver_->Step(block_unk, next_time, current_time_step);

  //// Updating vect_unk
  for (size_t i = 0; i < unk_size; i++) {
    auto &unk_i = *(vect_unk[i]);
    const mfem::Vector &bb = block_unk.GetBlock(i);
    unk_i = bb;
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template <class T, int DIM, class NLFI, class LHS_NLFI>
void TransientOperatorBase<T, DIM, NLFI, LHS_NLFI>::get_mass_coefficient(const mfem::Vector &u) {
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
template <class T, int DIM, class NLFI, class LHS_NLFI>
void TransientOperatorBase<T, DIM, NLFI, LHS_NLFI>::build_mass_matrix(
    const std::vector<mfem::Vector> &u_vect) {
  this->M_solver_.clear();
  for (int i = 0; i < u_vect.size(); i++) {
    if (M != nullptr) {
      delete M;
    }
    ////////////////
    // Mass matrix (constant)
    ////////////////
    M = new mfem::ParBilinearForm(this->fes_[i]);

    auto mass_coefficient = mfem::ConstantCoefficient(1.0);
    if (!this->isExplicit_) {
      M->AddDomainIntegrator(new mfem::MassIntegrator(mass_coefficient));
    } else {
      M->AddDomainIntegrator(
          new mfem::LumpedIntegrator(new mfem::MassIntegrator(mass_coefficient)));
    }
    M->Assemble(0);
    M->Finalize(0);

    Mmat = M->ParallelAssemble();
    std::unique_ptr<mfem::HypreParMatrix> Me(Mmat->EliminateRowsCols(this->ess_tdof_list_[i]));

    this->mass_matrix_solver_ =
        std::make_shared<LSolver>(this->mass_solver_, this->mass_solver_params_,
                                  this->mass_precond_, this->mass_precond_params_, *Mmat);
    this->M_solver_.emplace_back(this->mass_matrix_solver_->get_solver());
  }
}

/**
 * @brief Build the nonlinear form associated with LHS of the PDEs
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param dt
 * @param u
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void TransientOperatorBase<T, DIM, NLFI, LHS_NLFI>::build_lhs_nonlinear_form(
    const double dt, const std::vector<mfem::Vector> &u_vect) {
  if (LHS != nullptr) {
    delete LHS;
  }

  LHS = new mfem::ParBlockNonlinearForm(this->fes_);

  this->lhs_nlfi_ptr_ = set_lhs_nlfi_ptr(dt, u_vect);

  LHS->AddDomainIntegrator(this->lhs_nlfi_ptr_);
}

/**
 * @brief  Set the LHS NonLinearFormIntegrator
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param dt
 * @param u_vect
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
LHS_NLFI *TransientOperatorBase<T, DIM, NLFI, LHS_NLFI>::set_lhs_nlfi_ptr(
    const double dt, const std::vector<mfem::Vector> &u) {
  Catch_Time_Section("TransientOperatorBase::set_lhs_nlfi_ptr");

  std::vector<mfem::ParGridFunction> vun;
  for (int i = 0; i < u.size(); i++) {
    mfem::ParGridFunction un(this->fes_[i]);

    un.SetFromTrueDofs(u[i]);

    vun.emplace_back(un);
  }
  const Parameters &all_params = this->params_ - this->default_p_;

  LHS_NLFI *lhs_nlfi_ptr = new LHS_NLFI(vun, all_params, this->auxvariables_);

  return lhs_nlfi_ptr;
}

/**
 * @brief Set current dt, unk values - needed to compute action and Jacobian.
 *solution_coef
 * @param dt time-step
 * @param u unknown vector
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void TransientOperatorBase<T, DIM, NLFI, LHS_NLFI>::SetTransientParameters(
    const double dt, const std::vector<mfem::Vector> &u_vect) {
  Catch_Time_Section("TransientOperatorBase::SetTransientParameters");

  ////////////////////////////////////////////
  // Build the LHS of the PDEs
  ////////////////////////////////////////////
  this->build_lhs_nonlinear_form(dt, u_vect);

  ////////////////////////////////////////////
  //  Build the RHS of the PDEs
  ////////////////////////////////////////////
  this->build_rhs_nonlinear_form(dt, u_vect);

  ////////////////////////////////////////////
  // Build Newton Linear system
  ////////////////////////////////////////////
  if (reduced_oper != nullptr) {
    delete reduced_oper;
  }
  reduced_oper = new PhaseFieldReducedOperator(this->LHS, this->RHS, this->ess_tdof_list_);
  ////////////////////////////////////////////
  // Newton Solver
  ////////////////////////////////////////////
  this->SetNewtonAlgorithm(reduced_oper);
}
/**
 * @brief Compute the mass matrix and the non linear form with the solution at the previous time
 * step
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void TransientOperatorBase<T, DIM, NLFI, LHS_NLFI>::SetExplicitTransientParameters(
    const std::vector<mfem::Vector> &un_vect) {
  Catch_Time_Section("TransientOperatorBase::SetExplicitTransientParameters");
  ////////////////////////////////////////////
  // Variable mass matrix
  ////////////////////////////////////////////
  this->build_mass_matrix(un_vect);

  ////////////////////////////////////////////
  // PhaseField non linear form
  ////////////////////////////////////////////
  this->build_rhs_nonlinear_form(0., un_vect);
}

/**
 * @brief Set current dt, unk values - needed to compute action and Jacobian.
 *
 * @param dt time-step
 * @param u unknown vector
 * @param ess_tdof_list array of dofs
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void TransientOperatorBase<T, DIM, NLFI, LHS_NLFI>::SetConstantParameters(
    const double dt, const std::vector<mfem::Vector> &u_vect) {
  Catch_Time_Section("TransientOperatorBase::SetConstantParameters");
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
template <class T, int DIM, class NLFI, class LHS_NLFI>
void TransientOperatorBase<T, DIM, NLFI, LHS_NLFI>::Mult(const mfem::Vector &u,
                                                         mfem::Vector &du_dt) const {
  Catch_Time_Section("TransientOperatorBase::Mult");

  const auto sc = this->height_;
  mfem::Vector v(u.GetData(), sc);
  mfem::Vector dv_dt(du_dt.GetData(), sc);

  // Todo(cci) change with BlockVector
  std::vector<mfem::Vector> v_vect;
  v_vect.emplace_back(v);

  // Todo(cci) : try to do different because of not satisfying
  const_cast<TransientOperatorBase<T, DIM, NLFI, LHS_NLFI> *>(this)->SetExplicitTransientParameters(
      v_vect);

  // const mfem::Array<int> offsets = this->RHS->GetBlockOffsets();
  const int fes_size = this->block_trueOffsets_.Size() - 1;
  mfem::BlockVector bb(this->block_trueOffsets_);

  this->RHS->Mult(v, bb);

  // Source term
  mfem::BlockVector source_term(this->block_trueOffsets_);
  source_term = 0.0;
  if (!this->src_func_.empty()) {
    for (int i = 0; i < fes_size; ++i) {
      if (this->src_func_[i] != nullptr) {
        mfem::ParLinearForm *RHS = new mfem::ParLinearForm(this->fes_[i]);
        mfem::Vector &src_i = source_term.GetBlock(i);
        this->get_source_term(i, this->src_func_[i], src_i, RHS);

        delete RHS;
      }
    }
    bb -= source_term;
  }
  bb.Neg();

  mfem::BlockVector sol(this->block_trueOffsets_);
  sol = 0.;

  for (int i = 0; i < fes_size; ++i) {
    mfem::Vector &soli = sol.GetBlock(i);
    mfem::Vector &bb_i = bb.GetBlock(i);
    std::visit(
        [&](auto &&arg) {
          using TT = std::decay_t<decltype(arg)>;
          if constexpr (!std::is_same_v<TT, std::shared_ptr<std::monostate>>) {
            if (bb_i.Size() != arg->NumCols()) {
              std::string msg = "Size mismatch: bb_i.Size() [" + std::to_string(bb_i.Size()) +
                                "] != arg->NumCols() [" + std::to_string(arg->NumCols()) + "]";
              MFEM_ABORT(msg.c_str());
            }
            arg->Mult(bb_i, soli);
          }
        },
        this->M_solver_[i]);
  }
  dv_dt = sol;
}

/**
 * @brief  Solve the Backward-Euler equation: k = f(phi + dt*k, t), for the unknown k.
 *
 * @param dt current time step
 * @param u unknown vector
 * @param du_dt unkwon time derivative vector
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void TransientOperatorBase<T, DIM, NLFI, LHS_NLFI>::ImplicitSolve(const double dt,
                                                                  const mfem::Vector &u,
                                                                  mfem::Vector &du_dt) {
  Catch_Time_Section("TransientOperatorBase::ImplicitSolve");

  const auto sc = this->height_;
  mfem::Vector v(u.GetData(), sc);
  mfem::Vector dv_dt(du_dt.GetData(), sc);
  const int fes_size = this->block_trueOffsets_.Size() - 1;

  {
    Catch_Time_Section("ImplicitSolve::SetTransientParams");
    std::vector<mfem::Vector> v_vect;
    auto sc_1 = 0;
    auto sc_2 = sc / fes_size;
    for (int i = 0; i < fes_size; ++i) {
      mfem::Vector v_i(u.GetData() + sc_1, sc_2);
      sc_1 += sc_2;
      v_vect.emplace_back(v_i);
    }
    this->SetTransientParameters(dt, v_vect);
  }
  // Apply BCs
  {
    Catch_Time_Section("ImplicitSolve::ApplyBCs");
    auto sc_1 = 0;
    auto sc_2 = sc / fes_size;
    for (int i = 0; i < fes_size; ++i) {
      mfem::Vector v_i(u.GetData() + sc_1, sc_2);
      this->bcs_[i]->SetBoundaryConditions(v_i);
      sc_1 += sc_2;
    }
    reduced_oper->SetParameters(dt, &v);
  }
  // Source term
  mfem::BlockVector source_term(this->block_trueOffsets_);
  source_term = 0.0;
  {
    Catch_Time_Section("ImplicitSolve::SourceTerm");
    if (!this->src_func_.empty()) {
      for (int i = 0; i < fes_size; ++i) {
        if (this->src_func_[i] != nullptr) {
          mfem::ParLinearForm *RHS = new mfem::ParLinearForm(this->fes_[i]);
          mfem::Vector &src_i = source_term.GetBlock(i);
          this->get_source_term(i, this->src_func_[i], src_i, RHS);
          delete RHS;
        }
      }
    }
  }

  // source_term.Print();
  {
    Catch_Time_Section("ImplicitSolve::CallMult");
    this->newton_solver_->Mult(source_term, dv_dt);
    delete this->rhs_solver_;
  }

  MFEM_VERIFY(this->newton_solver_->GetConverged(), "Nonlinear solver did not converge.");
}

/**
 * @brief Overload the solver used to invert the mass matrix
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param SOLVER
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void TransientOperatorBase<T, DIM, NLFI, LHS_NLFI>::overload_mass_solver(VSolverType SOLVER) {
  this->mass_solver_ = SOLVER;
}

/**
 * @brief Overload the solver used to invert the mass matrix with its parameters
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param SOLVER
 * @param s_params
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void TransientOperatorBase<T, DIM, NLFI, LHS_NLFI>::overload_mass_solver(
    VSolverType SOLVER, const Parameters &s_params) {
  this->mass_solver_ = SOLVER;
  this->mass_solver_params_ = s_params;
}

/**
 * @brief Overload the preconditioner for the solver used to invert the mass matrix
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param PRECOND
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void TransientOperatorBase<T, DIM, NLFI, LHS_NLFI>::overload_mass_preconditioner(
    VSolverType PRECOND) {
  this->mass_precond_ = PRECOND;
}

/**
 * @brief Overload the preconditioner for the solver used to invert the mass matrix and its
 * parameters
 *
 * @tparam T
 * @tparam DIM
 * @tparam NLFI
 * @param PRECOND
 * @param p_params
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
void TransientOperatorBase<T, DIM, NLFI, LHS_NLFI>::overload_mass_preconditioner(
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
template <class T, int DIM, class NLFI, class LHS_NLFI>
void TransientOperatorBase<T, DIM, NLFI, LHS_NLFI>::set_default_mass_solver() {
  auto s_params = Parameters(Parameter("description", "Default Mass Solver"));
  auto p_params = Parameters(Parameter("description", "Default Mass preconditioner"));

  this->mass_solver_ = HypreSolverType::HYPRE_GMRES;
  this->mass_solver_params_ = s_params;
  this->mass_precond_ = HyprePreconditionerType::HYPRE_ILU;
  this->mass_precond_params_ = p_params;
}
/**
 * @brief Destroy the Phase Field Operator:: Phase Field Operator object
 *
 */
template <class T, int DIM, class NLFI, class LHS_NLFI>
TransientOperatorBase<T, DIM, NLFI, LHS_NLFI>::~TransientOperatorBase() {}

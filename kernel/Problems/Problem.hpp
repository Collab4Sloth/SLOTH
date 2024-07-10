/**
 * @file Problem.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Class used to defined a Problem objet
 * @version 0.1
 * @date 2024-05-25
 *
 * Copyright CEA (c) 2024
 *
 */

#pragma once
#include <memory>
#include <string>
#include <tuple>

#include "Convergence/PhysicalConvergence.hpp"
#include "Parameters/Parameter.hpp"
#include "PostProcessing/postprocessing.hpp"
#include "Variables/Variable.hpp"
#include "mfem.hpp" // NOLINT [no include the directory when naming mfem include file]

template <class OPE, class VAR, class PST>
class Problem {
 private:
  std::string name_{"UNKNOWN PROBLEM"};
  mfem::ODESolver* ode_solver_;
  OPE oper_;
  VAR variables_;
  PST& pst_;
  mfem::Vector unknown_;
  PhysicalConvergence convergence_;
  std::shared_ptr<std::function<double(const mfem::Vector&, double)> > analytical_solution_{
      nullptr};
  void set_ODE_solver(const TimeScheme::value& ode_solver);
  void do_time_step(mfem::Vector& unk, double& current_time, double current_time_step);
  std::tuple<bool, double> check_convergence(const mfem::Vector& unk, const mfem::Vector& prev_unk);

 public:
  Problem(const std::string& name, const OPE& oper, const VAR& variables, PST& pst,
          TimeScheme::value ode, const PhysicalConvergence& convergence, const Parameters& params);
  const std::string get_name();
  VAR get_problem_variables();
  void update();
  void initialize(const double& initial_time);
  std::tuple<bool, double, mfem::Vector> execute(const int& iter, double& current_time,
                                                 const double& current_time_step);
  void post_execute(const int& iter, const double& current_time_step, const double& current_time);
  void save(const int& iter, const double& current_time);
  void post_processing(const int& iter, const double& current_time_step,
                       const double& current_time);
  void finalize();
  ~Problem();
};

/**
 * @brief Construct a new Problem object
 *
 * @tparam OPE
 * @tparam SOLVER
 * @param oper
 * @param solver
 * @param convergence
 */
template <class OPE, class VAR, class PST>
Problem<OPE, VAR, PST>::Problem(const std::string& name, const OPE& oper, const VAR& variables,
                                PST& pst, TimeScheme::value ode,
                                const PhysicalConvergence& convergence, const Parameters& params)
    : name_(name), oper_(oper), variables_(variables), pst_(pst), convergence_(convergence) {
  this->set_ODE_solver(ode);
}

/**
 * @brief Set the ODE time marching
 *
 * @tparam OPE
 * @tparam VAR
 * @param ode_solver
 */
template <class OPE, class VAR, class PST>
void Problem<OPE, VAR, PST>::set_ODE_solver(const TimeScheme::value& ode_solver) {
  switch (ode_solver) {
    case TimeScheme::EulerExplicit: {
      this->ode_solver_ = new mfem::ForwardEulerSolver;
      break;
    }
    case TimeScheme::EulerImplicit: {
      this->ode_solver_ = new mfem::BackwardEulerSolver;
      break;
    }
    case TimeScheme::RungeKutta4: {
      this->ode_solver_ = new mfem::SDIRK33Solver;
      break;
    }
    default:
      throw std::runtime_error(
          "TimeDiscretization::set_ODE_solver: EulerImplicit, EulerExplicit, Rungekutta4 are "
          "available");
      break;
  }
}

/**
 * @brief Return the name of the problem
 *
 * @tparam OPE
 * @tparam VAR
 * @tparam PST
 * @return const std::string
 */
template <class OPE, class VAR, class PST>
const std::string Problem<OPE, VAR, PST>::get_name() {
  return this->name_;
}

/**
 * @brief Return the variables associated with the problem
 *
 * @tparam OPE
 * @tparam VAR
 * @tparam PST
 * @return VAR
 */
template <class OPE, class VAR, class PST>
VAR Problem<OPE, VAR, PST>::get_problem_variables() {
  return this->variables_;
}

/**
 * @brief Update the variables associated with the problem
 *
 * @tparam OPE
 * @tparam VAR
 * @tparam PST
 */
template <class OPE, class VAR, class PST>
void Problem<OPE, VAR, PST>::update() {
  // auto& var = this->variables_.get_variable("phi");
  auto& var = this->variables_.getIVariable(0);
  var.update(this->unknown_);
}

/**
 * @brief  Initialize the calculation : operator + ODE
 *
 * @tparam OPE
 * @tparam VAR
 * @tparam PST
 * @param initial_time
 */
template <class OPE, class VAR, class PST>
void Problem<OPE, VAR, PST>::initialize(const double& initial_time) {
  this->oper_.initialize(initial_time);
  // Call before the first call to step() or when the time-step is restarted
  this->ode_solver_->Init(this->oper_);
}

/**
 * @brief Run a time-step : calculation + check of convergence
 *
 * @tparam OPE
 * @tparam VAR
 * @tparam PST
 * @param iter
 * @param current_time
 * @param current_time_step
 * @return std::tuple<bool, double, mfem::Vector>
 */
template <class OPE, class VAR, class PST>
std::tuple<bool, double, mfem::Vector> Problem<OPE, VAR, PST>::execute(
    const int& iter, double& current_time, const double& current_time_step) {
  SlothInfo::print("   ============================== ");
  SlothInfo::print("   ==== Problem : ", this->name_);
  SlothInfo::print("   ============================== ");
  auto& var = this->variables_.getIVariable(0);
  auto unk = var.get_unknown();
  this->do_time_step(unk, current_time, current_time_step);

  bool is_converged = true;
  auto criterion = 0.;
  if (iter > 1) {
    const auto& prev_unk = var.get_last();
    const auto& [is_converged_iter, criterion_iter] = this->check_convergence(unk, prev_unk);
    is_converged = is_converged_iter;
    criterion = criterion_iter;
  }
  return std::make_tuple(is_converged, criterion, unk);
}

/**
 * @brief Compute Error and Energies associated with the problem if required
 *
 * @tparam OPE
 * @tparam VAR
 * @param iter
 * @param current_time_step
 * @param current_time
 * @param unk
 */
template <class OPE, class VAR, class PST>
void Problem<OPE, VAR, PST>::post_execute(const int& iter, const double& current_time_step,
                                          const double& current_time) {
  auto vv = this->variables_.getIVariable(0);
  auto unk = vv.get_unknown();
  auto solution = vv.get_analytical_solution();
  if (solution != nullptr) {
    auto solution_func = solution.get();

    this->oper_.ComputeError(iter, current_time_step, current_time, unk, *solution_func);
  }
  this->oper_.ComputeEnergies(iter, current_time_step, current_time, unk);
}

/**
 * @brief Fill CSV files, if required
 *
 * @tparam OPE
 * @tparam VAR
 */
template <class OPE, class VAR, class PST>
void Problem<OPE, VAR, PST>::finalize() {
  this->pst_.save_specialized(this->oper_.get_time_specialized());
}

/**
 * @brief Call the post_execute method of the given problem and saves variables according with PST
 *
 * @tparam OPE
 * @tparam VAR
 * @tparam PST
 * @param iter
 * @param current_time
 * @param current_time_step
 */
template <class OPE, class VAR, class PST>
void Problem<OPE, VAR, PST>::post_processing(const int& iter, const double& current_time,
                                             const double& current_time_step) {
  this->post_execute(iter, current_time_step, current_time);
  // Save for visualization
  this->save(iter, current_time);
}

/**
 * @brief Do a time-step by calling Step method of the ODE
 *
 * @tparam OPE
 * @tparam VAR
 * @param unk
 * @param current_time
 * @param current_time_step
 */
template <class OPE, class VAR, class PST>
void Problem<OPE, VAR, PST>::do_time_step(mfem::Vector& unk, double& current_time,
                                          double current_time_step) {
  this->ode_solver_->Step(unk, current_time, current_time_step);

  // Store the solution into a temporary mfem::Vector that will be used during updating stage, if
  // calculation converges
  this->unknown_ = unk;
}

/**
 * @brief Check convergence at the current iteration
 *
 * @tparam OPE
 * @tparam VAR
 * @param unk
 * @param prev_unk
 * @return std::tuple<bool, double>
 */
template <class OPE, class VAR, class PST>
std::tuple<bool, double> Problem<OPE, VAR, PST>::check_convergence(const mfem::Vector& unk,
                                                                   const mfem::Vector& prev_unk) {
  return this->convergence_.getPhysicalConvergence(unk, prev_unk);
}

template <class OPE, class VAR, class PST>
void Problem<OPE, VAR, PST>::save(const int& iter, const double& current_time) {
  auto vars = this->get_problem_variables();
  this->pst_.save_variables(vars, iter, current_time);
}

/**
 * @brief Destroy the Problem object
 *
 * @tparam OPE
 * @tparam SOLVER
 */
template <class OPE, class VAR, class PST>
Problem<OPE, VAR, PST>::~Problem() {}

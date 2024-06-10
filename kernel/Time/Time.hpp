/*
 * Copyright © CEA 2023
 *
 * \brief TIme discretization used by phase-field models
 *
 * \file TimeDiscretization.hpp
 * \author ci230846
 * \date 26/03/2023
 */

#pragma once
#include <functional>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include "PostProcessing/postprocessing.hpp"
#include "Utils/PhaseFieldOptions.hpp"
#include "Variables/Variable.hpp"
#include "mfem.hpp"


template <class PST, class OPE, class VAR>
class TimeDiscretization {
 private:
  mfem::ODESolver* ode_solver_;
  OPE oper_;
  VAR variables_;
  PST& pst_;
  double initial_time_{0.};
  double final_time_;
  double time_step_;
  bool compute_error_{false};
  bool compute_energies_{false};
  bool last_step_{false};
  std::shared_ptr<std::function<double(const mfem::Vector&, double)> > analytical_solution_{
      nullptr};

  void set_parameters(const Parameters& params);
  void set_ODE_solver(const std::string& ode_solver);
  void set_options(const Parameters& params);
  void initialize();
  double compute_real_time_step(const double& current_time);
  void print_time_step_info(const int& iter, const double& time, const double& dt_requested,
                            const double& dt_real);

 public:
  TimeDiscretization(const std::string& ode_solver, const OPE& oper, const Parameters& params,
                     const VAR& variables, PST& pst);
  void execute();
  ~TimeDiscretization();
};

/**
 * @brief Construct a new Time Discretization:: Time Discretization object
 *
 * @param ode_solver
 * @param unknown
 * @param with_save
 */
template <class PST, class OPE, class VAR>
TimeDiscretization<PST, OPE, VAR>::TimeDiscretization(const std::string& ode_solver,
                                                      const OPE& oper, const Parameters& params,
                                                      const VAR& variables, PST& pst)
    : oper_(oper), variables_(variables), pst_(pst) {
  this->set_parameters(params);
  this->set_ODE_solver(ode_solver);
  this->set_options(params);
}

/**
 * @brief
 *
 * @tparam T
 * @tparam DC
 * @tparam DIM
 * @tparam NLFI
 * @param ode_solver
 */
template <class PST, class OPE, class VAR>
void TimeDiscretization<PST, OPE, VAR>::set_ODE_solver(const std::string& ode_solver) {
  switch (TimeScheme::from(ode_solver)) {
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
 * @brief Set options of calculation
 *
 * @tparam PST
 * @tparam OPE
 * @tparam VAR
 * @param params
 */
template <class PST, class OPE, class VAR>
void TimeDiscretization<PST, OPE, VAR>::set_parameters(const Parameters& params) {
  this->initial_time_ = params.get_parameter_value("initial_time");
  this->final_time_ = params.get_parameter_value("final_time");
  this->time_step_ = params.get_parameter_value("time_step");
}

/**
 * @brief Set options of calculation, create CSV files if required
 *
 * @tparam PST
 * @tparam OPE
 * @tparam VAR
 * @param params
 */
template <class PST, class OPE, class VAR>
void TimeDiscretization<PST, OPE, VAR>::set_options(const Parameters& params) {
  this->compute_error_ = params.get_option_value("compute_error");
  if (this->compute_error_) {
    auto vv = variables_.get_variable("phi");
    this->analytical_solution_ = vv.get_analytical_solution();
    if (this->analytical_solution_ != nullptr) {
      const std::string& file = "error.csv";
      const std::string& headers = "Iter[-] Dt[s] Time[s] L2-error[-]";
      this->pst_.create_csv(file, headers);
    } else {
      throw std::runtime_error(
          "PhaseFieldOperatorBase<T, DIM, NLFI>::ComputeError() : analytical solution doesn't "
          "exist. "
          "Please check input data.");
    }
  }
  this->compute_energies_ = params.get_option_value("compute_energies");
  if (this->compute_energies_) {
    const std::string& file = "energy_density.csv";
    const std::string& headers = "Iter[-]  Dt[s] Time[s] Density[J.m-3]";
    this->pst_.create_csv(file, headers);
    const std::string& file1 = "energy_interface.csv";
    const std::string& headers1 = "Iter[-]  Dt[s] Time[s] Sigma[J.m-3]";
    this->pst_.create_csv(file1, headers1);
  }
}

/**
 * @brief Initialization stage
 *
 */
template <class PST, class OPE, class VAR>
void TimeDiscretization<PST, OPE, VAR>::initialize() {
  Timers timer_initialize("TimeDiscretization::initialize");
  UtilsForOutput::getInstance().get_indentation();
  timer_initialize.start();

  // this->oper_.SetTime(this->initial_time_);
  this->oper_.initialize(this->initial_time_);
  // Call before the first call to step() or when the time-step is restarted
  this->ode_solver_->Init(this->oper_);
  // Save at initial time
  this->pst_.save_variables(this->variables_, 0, this->initial_time_);

  timer_initialize.stop();
  UtilsForOutput::getInstance().update_timer("TimeDiscretization::initialize",timer_initialize);
  UtilsForOutput::getInstance().get_lose_indentation("TimeDiscretization::initialize");
}

/**
 * @brief Compute the current value of the time-step
 *
 * @return double
 */
template <class PST, class OPE, class VAR>
double TimeDiscretization<PST, OPE, VAR>::compute_real_time_step(const double& current_time) {
  auto dt_real = this->time_step_;
  if ((current_time + this->time_step_) + 0.5 * this->time_step_ > this->final_time_) {
    this->last_step_ = true;
    dt_real = this->final_time_ - current_time;
  }
  
  return dt_real;

  
}

/**
 * @brief Run the calculation
 *
 */
template <class PST, class OPE, class VAR>
void TimeDiscretization<PST, OPE, VAR>::execute() {
  //------Start profiling-------------------------;
  Timers timer_execute("TimeDiscretization::execute");
  UtilsForOutput::getInstance().get_indentation();
  timer_execute.start();
  //------------------------------------------------

  // Initialization
  this->initialize();
  // Variable& var2 = variables_.getIVariable(0);
  // TODO(cci) : passe par la map
  auto& var = variables_.get_variable("phi");
  // Not necessary : deja dans Operator au constructeur
  // this->oper_.initialize(var);
  auto current_time = this->initial_time_;
  // TODO(cci) : construire un block vector

  auto unk = var.get_unknown();

  //=============================================
  //           loop over time-step
  //=============================================


  for (auto iter = 1; !this->last_step_; iter++) {
    //------------
    // Time-step
    //------------
    auto current_time_step = compute_real_time_step(current_time);
    this->print_time_step_info(iter, current_time, this->time_step_, current_time_step);
    //------------
    // Solve
    //------------
    // TODO(cci) : passer le block vector

    this->ode_solver_->Step(unk, current_time, current_time_step);

    //-----------------
    // Update solution
    //-----------------
    // TODO(cci) : vers l'update avec une loop su rle block vector

    var.update(unk);

    //---------------------
    // Extra calculation
    //---------------------
    // oper.SetParameters(unk);

    
    if (this->compute_energies_) {
      this->oper_.ComputeEnergies(std::make_tuple(iter, current_time_step, current_time), unk);
    }
    
    

    if (this->compute_error_) {

      auto solution_func = this->analytical_solution_.get();

      this->oper_.ComputeError(std::make_tuple(iter, current_time_step, current_time), unk,
                               *solution_func);
      
    }
    
    //-------------------
    // Visualization
    //-------------------
    if (this->pst_.need_to_be_saved(iter)) {
      this->pst_.save_variables(this->variables_, iter, current_time);
    }
  }
  

  if (this->compute_error_) {
    this->pst_.export_csv("error.csv", this->oper_.get_l2_error());
  }

  
  
  if (this->compute_energies_) {
    this->pst_.export_csv("energy_density.csv", this->oper_.get_energy_density());
    this->pst_.export_csv("energy_interface.csv", this->oper_.get_energy_interface());
  }
  
  //=============================================
  //          end of loop over time-step
  //=============================================

  timer_execute.stop();
  UtilsForOutput::getInstance().update_timer("TimeDiscretization::execute",timer_execute);
  UtilsForOutput::getInstance().get_lose_indentation("TimeDiscretization::execute");
  
}

/**
 * @brief Print time-step information
 *
 * @tparam PST
 * @tparam OPE
 * @tparam VAR
 * @param iter
 * @param dt_requested
 * @param dt_real
 */
template <class PST, class OPE, class VAR>
void TimeDiscretization<PST, OPE, VAR>::print_time_step_info(const int& iter, const double& time,
                                                             const double& dt_requested,
                                                             const double& dt_real) {
  std::cout << " ============================== " << '\n';
  std::cout << " ==== Iteration : " << iter << '\n';
  std::cout << "      - time : " << time << '\n';
  std::cout << "      - time-step requested : " << dt_requested << '\n';
  std::cout << "      - time-step accepted  : " << dt_real << '\n';
  std::cout << " ============================== " << '\n';
}

/**
 * @brief Destroy the Time Discretization:: Time Discretization object
 *
 */
template <class PST, class OPE, class VAR>
TimeDiscretization<PST, OPE, VAR>::~TimeDiscretization() {}

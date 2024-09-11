/**
 * @file ProblemBase.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Base class for building Problems
 * @version 0.1
 * @date 2024-05-25
 *
 * Copyright CEA (c) 2024
 *
 */

#pragma once
#include <memory>
#include <optional>
#include <string>
#include <tuple>

#include "Convergence/PhysicalConvergence.hpp"
#include "Parameters/Parameter.hpp"
#include "PostProcessing/postprocessing.hpp"
#include "Variables/Variable.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

template <class VAR, class PST>
class ProblemBase {
 private:
  std::tuple<bool, double> check_convergence(const mfem::Vector& unk, const mfem::Vector& prev_unk);

 protected:
  std::string description_{"UNKNOWN PROBLEM"};
  VAR& variables_;
  VAR* auxvariables_;
  PST& pst_;
  mfem::Vector unknown_;
  PhysicalConvergence convergence_;

 public:
  ProblemBase(const std::string& name, VAR& variables, PST& pst,
              const PhysicalConvergence& convergence);
  ProblemBase(const std::string& name, VAR& variables, VAR& auxvariables, PST& pst,
              const PhysicalConvergence& convergence);

  std::string get_description();
  VAR get_problem_variables();

  /////////////////////////////////////////////////////

  virtual void initialize(const double& initial_time);

  /////////////////////////////////////////////////////
  std::tuple<bool, double, mfem::Vector> execute(const int& iter, double& next_time,
                                                 const double& current_time,
                                                 const double& current_time_step);

  virtual void do_time_step(mfem::Vector& unk, double& next_time, const double& current_time,
                            double current_time_step, const int iterp) = 0;

  virtual void post_execute(const int& iter, const double& current_time,
                            const double& current_time_step);

  /////////////////////////////////////////////////////

  void update();

  /////////////////////////////////////////////////////
  virtual void post_processing(const int& iter, const double& current_time,
                               const double& current_time_step);

  void save(const int& iter, const double& current_time);
  /////////////////////////////////////////////////////

  virtual void finalize();

  virtual ~ProblemBase();
};

template <class VAR, class PST>
ProblemBase<VAR, PST>::ProblemBase(const std::string& name, VAR& variables, PST& pst,
                                   const PhysicalConvergence& convergence)
    : description_(name),
      variables_(variables),
      auxvariables_(nullptr),
      pst_(pst),
      convergence_(convergence) {}

/**
 * @brief Construct a new ProblemBase<VAR, PST>::ProblemBase object
 *
 * @tparam VAR
 * @tparam PST
 * @param name
 * @param variables
 * @param auxvariables
 * @param pst
 * @param convergence
 */
template <class VAR, class PST>
ProblemBase<VAR, PST>::ProblemBase(const std::string& name, VAR& variables, VAR& auxvariables,
                                   PST& pst, const PhysicalConvergence& convergence)
    : description_(name),
      variables_(variables),
      auxvariables_(&auxvariables),
      pst_(pst),
      convergence_(convergence) {}

/**
 * @brief Return the name of the problem
 *
 * @tparam VAR
 * @tparam PST
 * @return const std::string
 */
template <class VAR, class PST>
std::string ProblemBase<VAR, PST>::get_description() {
  return this->description_;
}

/**
 * @brief  Return the variables associated with the problem
 *
 * @tparam VAR
 * @tparam PST
 * @return VAR
 */
template <class VAR, class PST>
VAR ProblemBase<VAR, PST>::get_problem_variables() {
  return this->variables_;
}

/**
 * @brief Action done after execute method
 *
 * @tparam VAR
 * @tparam PST
 * @param iter
 * @param current_time
 * @param current_time_step
 */
template <class VAR, class PST>
void ProblemBase<VAR, PST>::post_execute(const int& iter, const double& current_time,
                                         const double& current_time_step) {}

/**
 * @brief Update the variables associated with the problem
 *
 * @tparam VAR
 * @tparam PST
 */
template <class VAR, class PST>
void ProblemBase<VAR, PST>::update() {
  // auto& var = this->variables_.get_variable("phi");
  // TODO(cci) généraliser pour avoir plusieurs variables (eg CALPHAD)
  auto& var = this->variables_.getIVariable(0);
  var.update(this->unknown_);
}

/**
 * @brief Initialize the problem
 *
 * @tparam VAR
 * @tparam PST
 * @param initial_time
 */
template <class VAR, class PST>
void ProblemBase<VAR, PST>::initialize(const double& initial_time) {}

/**
 * @brief Run a time-step : calculation + check of convergence
 *
 * @tparam VAR
 * @tparam PST
 * @param iter
 * @param current_time
 * @param current_time_step
 * @return std::tuple<bool, double, mfem::Vector>
 */
template <class VAR, class PST>
std::tuple<bool, double, mfem::Vector> ProblemBase<VAR, PST>::execute(
    const int& iter, double& next_time, const double& current_time,
    const double& current_time_step) {
  int rank = mfem::Mpi::WorldRank();
  if (rank == 0) {
    SlothInfo::verbose("   ============================== ");
    SlothInfo::verbose("   ==== Problem : ", this->description_);
    SlothInfo::verbose("   ============================== ");
  }
  auto& var = this->variables_.getIVariable(0);
  auto unk = var.get_unknown();
  this->do_time_step(unk, next_time, current_time, current_time_step, iter);

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
 * @brief Finalize the problem
 *
 * @tparam VAR
 * @tparam PST
 */
template <class VAR, class PST>
void ProblemBase<VAR, PST>::finalize() {
  // finalize....
}

/**
 * @brief Save variables of the problem (see save method for other actions)
 *
 * @tparam VAR
 * @tparam PST
 * @param iter
 * @param current_time
 * @param current_time_step
 */
template <class VAR, class PST>
void ProblemBase<VAR, PST>::post_processing(const int& iter, const double& current_time,
                                            const double& current_time_step) {
  // Save for visualization
  this->save(iter, current_time);
}

/**
 * @brief Main actions done during a time-step
 *
 * @tparam VAR
 * @tparam PST
 * @param unk
 * @param current_time
 * @param current_time_step
 */
template <class VAR, class PST>
void ProblemBase<VAR, PST>::do_time_step(mfem::Vector& unk, double& next_time,
                                         const double& current_time, double current_time_step,
                                         const int iter) {}

/**
 * @brief Check convergence at the current iteration
 *
 * @tparam VAR
 * @tparam PST
 * @param unk
 * @param prev_unk
 * @return std::tuple<bool, double>
 */
template <class VAR, class PST>
std::tuple<bool, double> ProblemBase<VAR, PST>::check_convergence(const mfem::Vector& unk,
                                                                  const mfem::Vector& prev_unk) {
  return this->convergence_.getPhysicalConvergence(unk, prev_unk);
}

/**
 * @brief Save variables
 *
 * @tparam VAR
 * @tparam PST
 * @param iter
 * @param current_time
 */
template <class VAR, class PST>
void ProblemBase<VAR, PST>::save(const int& iter, const double& current_time) {
  auto vars = this->get_problem_variables();
  this->pst_.save_variables(vars, iter, current_time);
}

/**
 * @brief Destroy the ProblemBase<VAR, PST>::ProblemBase object
 *
 * @tparam VAR
 * @tparam PST
 */
template <class VAR, class PST>
ProblemBase<VAR, PST>::~ProblemBase() {}

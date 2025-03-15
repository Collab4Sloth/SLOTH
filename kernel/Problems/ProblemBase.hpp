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
#include <functional>
#include <list>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

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
  std::vector<VAR*> auxvariables_;
  PST& pst_;
  const std::list<int> pop_elem_;
  std::vector<mfem::Vector> unknown_;
  PhysicalConvergence convergence_;

 public:

  template <class... Args>
  ProblemBase(const std::string& name, VAR& variables, PST& pst,
              const PhysicalConvergence& convergence, std::list<int> pop_elem,
              Args&&... auxvariables);

  template <class... Args>
  ProblemBase(const std::string& name, VAR& variables, PST& pst,
              const PhysicalConvergence& convergence, Args&&... auxvariables);

  std::string get_description();
  VAR get_problem_variables();

  /////////////////////////////////////////////////////

  virtual void initialize(const double& initial_time);

  /////////////////////////////////////////////////////
  std::tuple<bool, double, std::vector<mfem::Vector>> execute(const int& iter, double& next_time,
                                                              const double& current_time,
                                                              const double& current_time_step);

  virtual void do_time_step(double& next_time, const double& current_time, double current_time_step,
                            const int iter, std::vector<std::unique_ptr<mfem::Vector>>& unks,
                            const std::vector<std::vector<std::string>>& unks_info) = 0;

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

/**
 * @brief Construct a new Problem Base< V A R,  P S T>:: Problem Base object
 *
 * @tparam VAR
 * @tparam PST
 * @tparam Args
 * @param name
 * @param variables
 * @param pst
 * @param convergence
 * @param pop_elem
 * @param auxvariables
 */
template <class VAR, class PST>
template <class... Args>
ProblemBase<VAR, PST>::ProblemBase(const std::string& name, VAR& variables, PST& pst,
                                   const PhysicalConvergence& convergence, std::list<int> pop_elem,
                                   Args&&... auxvariables)
    : description_(name),
      variables_(variables),
      pst_(pst),
      pop_elem_(pop_elem),
      convergence_(convergence) {
  if constexpr (sizeof...(auxvariables) == 0) {
    this->auxvariables_.resize(0);
  } else {
    this->auxvariables_ = {&auxvariables...};

    if (pop_elem.size() > 0) {
      for (const auto& pp : pop_elem_) {
        if (pp < this->auxvariables_.size()) {
          this->auxvariables_.erase(std::next(this->auxvariables_.begin(), pp),
                                    std::next(this->auxvariables_.begin(), pp + 1));
        } else {
          std::string msg = "ProblemBase: Error with index use to pop element ";
          mfem::mfem_error(msg.c_str());
        }
      }
    }
  }
}

/**
 * @brief Construct a new Problem Base< V A R,  P S T>:: Problem Base object
 *
 * @tparam VAR
 * @tparam PST
 * @tparam Args
 * @param name
 * @param variables
 * @param pst
 * @param convergence
 * @param auxvariables
 */
template <class VAR, class PST>
template <class... Args>
ProblemBase<VAR, PST>::ProblemBase(const std::string& name, VAR& variables, PST& pst,
                                   const PhysicalConvergence& convergence, Args&&... auxvariables)
    : description_(name), variables_(variables), pst_(pst), convergence_(convergence) {
  if constexpr (sizeof...(auxvariables) == 0) {
    this->auxvariables_.resize(0);
  } else {
    this->auxvariables_ = {&auxvariables...};
  }
}

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
  MFEM_VERIFY(this->variables_.get_variables_number() == this->unknown_.size(),
              "Error while updating the variables: the number of variables must be equal to the "
              "size of the unknown vector");
  for (std::size_t i = 0; i < this->variables_.get_variables_number(); i++) {
    auto& var = this->variables_.getIVariable(i);
    var.update(this->unknown_[i]);
  }
  this->unknown_.clear();
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
std::tuple<bool, double, std::vector<mfem::Vector>> ProblemBase<VAR, PST>::execute(
    const int& iter, double& next_time, const double& current_time,
    const double& current_time_step) {
  int rank = mfem::Mpi::WorldRank();
  if (rank == 0) {
    SlothInfo::verbose("   ============================== ");
    SlothInfo::verbose("   ==== Problem : ", this->description_);
    SlothInfo::verbose("   ============================== ");
  }

  std::vector<std::unique_ptr<mfem::Vector>> vect_unk;
  std::vector<std::vector<std::string>> vect_unk_info;

  size_t num_vars = this->variables_.getVariables().size();
  vect_unk.reserve(num_vars);
  vect_unk_info.reserve(num_vars);

  for (auto& var : this->variables_.getVariables()) {
    vect_unk.push_back(std::make_unique<mfem::Vector>(var.get_unknown()));

    auto& unk_info = vect_unk_info.emplace_back(var.get_additional_variable_info());
    unk_info.insert(unk_info.begin(), var.getVariableName());
  }

  this->do_time_step(next_time, current_time, current_time_step, iter, vect_unk, vect_unk_info);

  bool is_converged = true;
  auto criterion = 0.;
  std::vector<mfem::Vector> vunk;
  if (iter > 1) {
    int j = 0;
    for (auto& var : this->variables_.getVariables()) {
      const auto& prev_unk = var.get_last();
      mfem::Vector& unk = *(vect_unk[j]);
      vunk.emplace_back(unk);
      const auto& [is_converged_iter, criterion_iter] = this->check_convergence(unk, prev_unk);
      is_converged = is_converged_iter;
      criterion = criterion_iter;
    }
  }
  return std::make_tuple(is_converged, criterion, vunk);
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
void ProblemBase<VAR, PST>::do_time_step(double& next_time, const double& current_time,
                                         double current_time_step, const int iter,
                                         std::vector<std::unique_ptr<mfem::Vector>>& unks,
                                         const std::vector<std::vector<std::string>>& unks_info) {}

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

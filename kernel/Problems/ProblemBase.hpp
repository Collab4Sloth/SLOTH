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
  void check_convergence(const std::vector<std::unique_ptr<mfem::Vector>>& unks);

 protected:
  std::string name_{"Unnamed problem"};
  VAR& variables_;
  std::vector<VAR*> auxvariables_;
  PST& pst_;
  const std::list<int> pop_elem_;
  std::vector<mfem::Vector> unknown_;
  std::shared_ptr<PhysicalConvergence> convergence_;
  std::vector<std::tuple<std::string, bool, double>> var_convergence_;

 public:
  template <class... Args>
  ProblemBase(const std::string& name, VAR& variables, PST& pst, PhysicalConvergence& convergence,
              std::list<int> pop_elem, Args&&... auxvariables);

  template <class... Args>
  ProblemBase(const std::string& name, VAR& variables, PST& pst, std::list<int> pop_elem,
              Args&&... auxvariables);

  template <class... Args>
  ProblemBase(const std::string& name, VAR& variables, PST& pst, PhysicalConvergence& convergence,
              Args&&... auxvariables);

  template <class... Args>
  ProblemBase(const std::string& name, VAR& variables, PST& pst, Args&&... auxvariables);

  std::string get_name();
  VAR get_problem_variables();

  std::vector<std::tuple<std::string, bool, double>> get_convergence();

  /////////////////////////////////////////////////////

  virtual void initialize(const double& initial_time);

  /////////////////////////////////////////////////////
  void execute(const int& iter, double& next_time, const double& current_time,
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
                                   PhysicalConvergence& convergence, std::list<int> pop_elem,
                                   Args&&... auxvariables)
    : name_(name), variables_(variables), pst_(pst), pop_elem_(pop_elem) {
  this->convergence_ = std::make_shared<PhysicalConvergence>(convergence);
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
 * @param pop_elem
 * @param auxvariables
 */
template <class VAR, class PST>
template <class... Args>
ProblemBase<VAR, PST>::ProblemBase(const std::string& name, VAR& variables, PST& pst,
                                   std::list<int> pop_elem, Args&&... auxvariables)
    : name_(name), variables_(variables), pst_(pst), pop_elem_(pop_elem), convergence_(nullptr) {
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
                                   PhysicalConvergence& convergence, Args&&... auxvariables)
    : name_(name), variables_(variables), pst_(pst) {
  this->convergence_ = std::make_shared<PhysicalConvergence>(convergence);
  if constexpr (sizeof...(auxvariables) == 0) {
    this->auxvariables_.resize(0);
  } else {
    this->auxvariables_ = {&auxvariables...};
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
 * @param auxvariables
 */
template <class VAR, class PST>
template <class... Args>
ProblemBase<VAR, PST>::ProblemBase(const std::string& name, VAR& variables, PST& pst,
                                   Args&&... auxvariables)
    : name_(name), variables_(variables), pst_(pst), convergence_(nullptr) {
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
std::string ProblemBase<VAR, PST>::get_name() {
  return this->name_;
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
void ProblemBase<VAR, PST>::execute(const int& iter, double& next_time, const double& current_time,
                                    const double& current_time_step) {
  int rank = mfem::Mpi::WorldRank();
  if (rank == 0) {
    SlothInfo::verbose("   ============================== ");
    SlothInfo::verbose("   ==== Problem : ", this->name_);
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

  if (this->convergence_ != nullptr && iter > 1) {
    this->check_convergence(vect_unk);
  }
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
 * @brief Check convergence of variables at the current time-step
 *
 * @tparam VAR
 * @tparam PST
 * @param unks
 */
template <class VAR, class PST>
void ProblemBase<VAR, PST>::check_convergence(
    const std::vector<std::unique_ptr<mfem::Vector>>& unks) {
  this->var_convergence_.clear();

  int j = -1;
  for (auto& var : this->variables_.getVariables()) {
    j++;
    const auto& prev_unk = var.get_last();
    mfem::Vector& unk = *(unks[j]);
    const auto& [has_converged, criterion] =
        this->convergence_->getPhysicalConvergence(unk, prev_unk);
    this->var_convergence_.emplace_back(
        std::make_tuple(var.getVariableName(), has_converged, criterion));
  }
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
 * @brief Return for each variable, a tuple (string,string, boolean, double) where the boolean is
 * the flag indicating if the criterion is reached or not, and the double is the criterion reached.
 * The first string is the name of the problem, and the second the name of the variable.
 *
 * @tparam VAR
 * @tparam PST
 * @return std::vector<std::tuple<std::string, bool, double>>
 */
template <class VAR, class PST>
std::vector<std::tuple<std::string, bool, double>> ProblemBase<VAR, PST>::get_convergence() {
  return this->var_convergence_;
}

/**
 * @brief Destroy the ProblemBase<VAR, PST>::ProblemBase object
 *
 * @tparam VAR
 * @tparam PST
 */
template <class VAR, class PST>
ProblemBase<VAR, PST>::~ProblemBase() {}

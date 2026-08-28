/**
 * @anchor ProblemBase
 * @file ProblemBase.tpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Base class for building Problems
 * @version 0.1
 * @date 2025-09-05
 *
 *
 * @copyright CEA (C) 2025
 *
 * This file is part of SLOTH.
 *
 * SLOTH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SLOTH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
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

#include "Coefficients/Coefficients.hpp"
#include "Convergence/PhysicalConvergence.hpp"
#include "Options/ProblemsOptions.hpp"
#include "Parameters/Parameter.hpp"
#include "PostProcessing/postprocessing.hpp"
#include "Utils/UtilsForData.hpp"
#include "Variables/Variable.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Constructs a new ProblemBase object.
 *
 * @tparam VAR Type representing the Variables.
 * @tparam PST Type representing the PostProcessing.
 * @tparam Args Types of auxiliary variables passed variadically.
 *
 * @param name Name of the problem.
 * @param variables Reference to the Variables object.
 * @param pst Reference to the PostProcessing object.
 * @param convergence Reference to the Convergence  object.
 * @param auxvariables Variadic pack of auxiliary variables.
 * @param Coeff Vector of coefficients.
 */
template <class VAR, class PST>
template <PbVar<VAR>... Args>
ProblemBase<VAR, PST>::ProblemBase(const std::string& name, VAR& variables,
                                   const std::vector<Coefficients>& Coeff, Convergence& convergence,
                                   PST& pst, Args&&... auxvariables)
    : name_(name), variables_(variables), coefficients_(Coeff), pst_(pst) {
  this->convergence_ = std::make_shared<Convergence>(convergence);
  if constexpr (sizeof...(auxvariables) == 0) {
    this->auxvariables_.clear();
  } else {
    this->auxvariables_ = {&auxvariables...};
  }

  if (this->has_pst()) {
    this->problem_post_processing_directory_ = this->get_pst().get_post_processing_directory();
  }
}
/**
 * @brief Constructs a new ProblemBase object.
 *
 * @tparam VAR Type representing the Variables.
 * @tparam PST Type representing the PostProcessing.
 * @tparam Args Types of auxiliary variables passed variadically.
 *
 * @param name Name of the problem.
 * @param variables Reference to the Variables object.
 * @param convergence Reference to the Convergence  object.
 * @param auxvariables Variadic pack of auxiliary variables.
 * @param Coeff Vector of coefficients.
 */
template <class VAR, class PST>
template <PbVar<VAR>... Args>
ProblemBase<VAR, PST>::ProblemBase(const std::string& name, VAR& variables,
                                   const std::vector<Coefficients>& Coeff, Convergence& convergence,
                                   Args&&... auxvariables)
    : name_(name), variables_(variables), coefficients_(Coeff), pst_(std::nullopt) {
  this->convergence_ = std::make_shared<Convergence>(convergence);
  if constexpr (sizeof...(auxvariables) == 0) {
    this->auxvariables_.clear();
  } else {
    this->auxvariables_ = {&auxvariables...};
  }
}

/**
 * @brief Constructs a new ProblemBase object.
 *
 * @tparam VAR Type representing the Variables.
 * @tparam PST Type representing the PostProcessing.
 * @tparam Args Types of auxiliary variables passed variadically.
 * @tparam CArgs Types of coefficients passed variadically.
 *
 * @param name Name of the problem.
 * @param variables Reference to the Variables object.
 * @param pst Reference to the PostProcessing object.
 * @param auxvariables Variadic pack of auxiliary variables.
 * @param Coeff Variadic pack of coefficients.
 */
template <class VAR, class PST>
template <PbVar<VAR>... Args>
ProblemBase<VAR, PST>::ProblemBase(const std::string& name, VAR& variables,
                                   const std::vector<Coefficients>& Coeff, PST& pst,
                                   Args&&... auxvariables)
    : name_(name), variables_(variables), coefficients_(Coeff), pst_(pst), convergence_(nullptr) {
  if constexpr (sizeof...(auxvariables) == 0) {
    this->auxvariables_.clear();
  } else {
    this->auxvariables_ = {&auxvariables...};
  }

  if (this->has_pst()) {
    this->problem_post_processing_directory_ = this->get_pst().get_post_processing_directory();
  }
}
/**
 * @brief Constructs a new ProblemBase object.
 *
 * @tparam VAR Type representing the Variables.
 * @tparam PST Type representing the PostProcessing.
 * @tparam Args Types of auxiliary variables passed variadically.
 * @tparam CArgs Types of coefficients passed variadically.
 *
 * @param name Name of the problem.
 * @param variables Reference to the Variables object.
 * @param auxvariables Variadic pack of auxiliary variables.
 * @param Coeff Variadic pack of coefficients.
 */
template <class VAR, class PST>
template <PbVar<VAR>... Args>
ProblemBase<VAR, PST>::ProblemBase(const std::string& name, VAR& variables,
                                   const std::vector<Coefficients>& Coeff, Args&&... auxvariables)
    : name_(name),
      variables_(variables),
      coefficients_(Coeff),
      pst_(std::nullopt),
      convergence_(nullptr) {
  if constexpr (sizeof...(auxvariables) == 0) {
    this->auxvariables_.clear();
  } else {
    this->auxvariables_ = {&auxvariables...};
  }
}

/**
 * @brief Constructs a new ProblemBase object.
 *
 * @tparam VAR Type representing the Variables.
 * @tparam PST Type representing the PostProcessing.
 * @tparam Args Types of auxiliary variables passed variadically.
 *
 * @param name Name of the problem.
 * @param variables Reference to the Variables object.
 * @param pst Reference to the PostProcessing object.
 * @param convergence Reference to the Convergence  object.
 * @param auxvariables Variadic pack of auxiliary variables.
 */
template <class VAR, class PST>
template <PbVar<VAR>... Args>
ProblemBase<VAR, PST>::ProblemBase(const std::string& name, VAR& variables,
                                   Convergence& convergence, PST& pst, Args&&... auxvariables)
    : name_(name), variables_(variables), pst_(pst) {
  this->convergence_ = std::make_shared<Convergence>(convergence);
  if constexpr (sizeof...(auxvariables) == 0) {
    this->auxvariables_.clear();
  } else {
    this->auxvariables_ = {&auxvariables...};
  }
  if (this->has_pst()) {
    this->problem_post_processing_directory_ = this->get_pst().get_post_processing_directory();
  }
}
/**
 * @brief Constructs a new ProblemBase object.
 *
 * @tparam VAR Type representing the Variables.
 * @tparam PST Type representing the PostProcessing.
 * @tparam Args Types of auxiliary variables passed variadically.
 *
 * @param name Name of the problem.
 * @param variables Reference to the Variables object.
 * @param convergence Reference to the Convergence  object.
 * @param auxvariables Variadic pack of auxiliary variables.
 */
template <class VAR, class PST>
template <PbVar<VAR>... Args>
ProblemBase<VAR, PST>::ProblemBase(const std::string& name, VAR& variables,
                                   Convergence& convergence, Args&&... auxvariables)
    : name_(name), variables_(variables), pst_(std::nullopt) {
  this->convergence_ = std::make_shared<Convergence>(convergence);
  if constexpr (sizeof...(auxvariables) == 0) {
    this->auxvariables_.clear();
  } else {
    this->auxvariables_ = {&auxvariables...};
  }
}

/**
 * @brief Constructs a new ProblemBase object.
 *
 * @tparam VAR Type representing the Variables.
 * @tparam PST Type representing the PostProcessing.
 * @tparam Args Types of auxiliary variables passed variadically.
 *
 * @param name Name of the problem.
 * @param variables Reference to the Variables object.
 * @param pst Reference to the PostProcessing object.
 * @param auxvariables Variadic pack of auxiliary variables.
 */
template <class VAR, class PST>
template <PbVar<VAR>... Args>
ProblemBase<VAR, PST>::ProblemBase(const std::string& name, VAR& variables, PST& pst,
                                   Args&&... auxvariables)
    : name_(name), variables_(variables), pst_(pst), convergence_(nullptr) {
  if constexpr (sizeof...(auxvariables) == 0) {
    this->auxvariables_.clear();
  } else {
    this->auxvariables_ = {&auxvariables...};
  }
  if (this->has_pst()) {
    this->problem_post_processing_directory_ = this->get_pst().get_post_processing_directory();
  }
}
/**
 * @brief Constructs a new ProblemBase object.
 *
 * @tparam VAR Type representing the Variables.
 * @tparam PST Type representing the PostProcessing.
 * @tparam Args Types of auxiliary variables passed variadically.
 *
 * @param name Name of the problem.
 * @param variables Reference to the Variables object.
 * @param auxvariables Variadic pack of auxiliary variables.
 */
template <class VAR, class PST>
template <PbVar<VAR>... Args>
ProblemBase<VAR, PST>::ProblemBase(const std::string& name, VAR& variables, Args&&... auxvariables)
    : name_(name), variables_(variables), pst_(std::nullopt), convergence_(nullptr) {
  if constexpr (sizeof...(auxvariables) == 0) {
    this->auxvariables_.clear();
  } else {
    this->auxvariables_ = {&auxvariables...};
  }
}

/**
 *
 * @brief Returns the name of the problem
 *
 * @tparam VAR Type representing the Variables.
 * @tparam PST Type representing the PostProcessing.
 * @return std::string The name of the problem.
 */
template <class VAR, class PST>
std::string ProblemBase<VAR, PST>::get_name() {
  return this->name_;
}

/**
 *
 * @brief Set the name of the problem
 *
 * @tparam VAR Type representing the Variables.
 * @tparam PST Type representing the PostProcessing.
 * @param std::string The name of the problem.
 */
template <class VAR, class PST>
void ProblemBase<VAR, PST>::set_name(const std::string& new_name) {
  this->name_ = new_name;
}

/**
 *
 * @brief Set the name for global specialized CSV file of the problem
 *
 * @tparam VAR Type representing the Variables.
 * @tparam PST Type representing the PostProcessing.
 * @param std::string The name of or global specialized CSV file of the problem
 */
template <class VAR, class PST>
void ProblemBase<VAR, PST>::set_name_save_specialized(const std::string& new_name) {
  this->name_save_specialized_ = new_name;
}

/**
 * @brief  Return the variables associated with the problem
 *
 * @tparam VAR Type representing the Variables.
 * @tparam PST Type representing the PostProcessing.
 * @return VAR& The list of Variables of the problem.
 */
template <class VAR, class PST>
VAR& ProblemBase<VAR, PST>::get_problem_variables() {
  return this->variables_;
}

/**
 * @brief Return the post-processing directory for the problem
 *
 * @tparam VAR Type representing the Variables.
 * @tparam PST Type representing the PostProcessing.
 * @return VAR The list of Variables of the problem.
 */
template <class VAR, class PST>
std::string ProblemBase<VAR, PST>::get_problem_post_processing_directory() {
  return this->problem_post_processing_directory_;
}

/**
 * @brief Update the variables associated with the problem
 *
 * @tparam VAR Type representing the problem Variables.
 * @tparam PST Type representing the post-processing.
 */
template <class VAR, class PST>
void ProblemBase<VAR, PST>::update() {
  MFEM_VERIFY(this->variables_.get_variables_number() == this->unknown_.size(),
              "Error while updating the variables: the number of variables must be equal to the "
              "size of the unknown vector");
  for (std::size_t i = 0; i < this->variables_.get_variables_number(); i++) {
    auto& var = this->variables_[i];
    var.update(this->unknown_[i]);
  }
  this->unknown_.clear();
}

/**
 * @brief Execute step
 *
 * @tparam VAR Type representing the problem Variables.
 * @tparam PST Type representing the post-processing.
 * @param iter Current iteration number.
 * @param next_time Next simulation time (Current time+ current time step).
 * @param current_time Current simulation time.
 * @param current_time_step Current time-step.
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
 * @tparam VAR Type representing the problem Variables.
 * @tparam PST Type representing the post-processing.
 */
template <class VAR, class PST>
void ProblemBase<VAR, PST>::finalize() {
  // finalize....
}

/**
 * @brief Save variables of the problem (see save method for other actions)
 *
 * @tparam VAR Type representing the problem Variables.
 * @tparam PST Type representing the post-processing.
 * @param iter Current iteration number.
 * @param current_time Current simulation time.
 * @param current_time_step Current time-step.
 */
template <class VAR, class PST>
void ProblemBase<VAR, PST>::post_processing(const int& iter, const double& current_time) {
  // Save for visualization
  this->save(iter, current_time);
}

/**
 * @brief Main actions done during a time-step
 *
 * @tparam VAR Type representing the problem Variables.
 * @tparam PST Type representing the post-processing.
 * @param next_time Next simulation time (Current time+ current time step).
 * @param current_time Current simulation time.
 * @param current_time_step Current time-step.
 * @param iter Current iteration number.
 * @param unks Unknown at the current time-step.
 * @param unks_info Additionnal informations associated with the unkwon of the Problem.
 */
// template <class VAR, class PST>
// void ProblemBase<VAR, PST>::do_time_step(double& next_time, const double& current_time,
//                                          double current_time_step, const int iter,
//                                          std::vector<std::unique_ptr<mfem::Vector>>& unks,
//                                          const std::vector<std::vector<std::string>>& unks_info)
//                                          {}

/**
 * @brief Check convergence of variables at the current time-step
 *
 * @tparam VAR Type representing the problem Variables.
 * @tparam PST Type representing the post-processing.
 * @param unks Unknown at the current time-step.
 */
template <class VAR, class PST>
void ProblemBase<VAR, PST>::check_convergence(
    const std::vector<std::unique_ptr<mfem::Vector>>& unks) {
  this->var_convergence_.clear();

  std::vector<PhysicalConvergence> vect_convergence = this->convergence_->getPhysicalConvergence();

  int j = -1;
  for (auto& var : this->variables_.getVariables()) {
    j++;
    const auto& prev_unk = var.get_last();
    mfem::Vector& unk = *(unks[j]);
    const auto& [local_converged, local_criterion] =
        vect_convergence[j].getPhysicalConvergence(unk, prev_unk);

    int local_conv_int = static_cast<int>(local_converged);
    int global_conv_int = 0;
    double global_criterion = 0.0;

    MPI_Allreduce(&local_conv_int, &global_conv_int, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    MPI_Allreduce(&local_criterion, &global_criterion, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    // int rank = mfem::Mpi::WorldRank();
    // if (rank == 0) {
    bool has_converged = static_cast<bool>(global_conv_int);
    this->var_convergence_.emplace_back(
        std::make_tuple(var.getVariableName(), has_converged, global_criterion));
    // }
  }
}

/**
 * @brief Save variables
 *
 * @tparam VAR Type representing the problem Variables.
 * @tparam PST Type representing the post-processing.
 * @param iter Current iteration number.
 * @param current_time Current simulation time.
 */
template <class VAR, class PST>
void ProblemBase<VAR, PST>::save(const int iter, const double& current_time) {
  if (this->has_pst()) {
    auto vars = this->get_problem_variables();
    this->get_pst().save_variables(vars, iter, current_time);
  }
}

/**
 * @brief Return for each variable, a tuple (string,string, boolean, double) where the boolean is
 * the flag indicating if the criterion is reached or not, and the double is the criterion reached.
 * The first string is the name of the problem, and the second the name of the variable.
 *
 * @tparam VAR Type representing the problem Variables.
 * @tparam PST Type representing the post-processing.
 * @return std::vector<std::tuple<std::string, bool, double>>
 */
template <class VAR, class PST>
std::vector<std::tuple<std::string, bool, double>> ProblemBase<VAR, PST>::get_convergence() {
  return this->var_convergence_;
}

/**
 * @brief Define the geometry of the problem (Cartesian, Axisymmetric)
 *
 * @tparam VAR Type representing the problem Variables.
 * @tparam PST Type representing the post-processing.
 * @param geometry
 */
template <class VAR, class PST>
void ProblemBase<VAR, PST>::setGeometry(Geometry geometry) {
  this->geometry_ = geometry;
}

/**
 * @brief Apply the pre-refinement loop on this Problem's initial
 *        condition, if an AMR strategy is attached.
 *
 * @details No-op if `set_amr()` was never called (`amr_ == nullptr`).
 *
 * @tparam VAR Type representing the problem Variables.
 * @tparam PST Type representing the post-processing.
 * @param all_vars Variables of every Problem in the enclosing Coupling,
 *                used to resynchronize the whole coupling after any mesh
 *                mutation (see `Coupling::CollectAllVariables()`).
 */
template <class VAR, class PST>
void ProblemBase<VAR, PST>::initialize_amr(const std::vector<VAR*>& all_vars) {
  if (this->amr_ != nullptr) {
    this->amr_->InitialRefine(this->variables_, all_vars);
  }
}

/**
 * @brief Apply one derefine-then-refine AMR pass for this Problem, if an
 *        AMR strategy is attached.
 *
 * @details No-op if `set_amr()` was never called (`amr_ == nullptr`).
 *
 * @tparam VAR Type representing the problem Variables.
 * @tparam PST Type representing the post-processing.
 * @param all_vars Variables of every Problem in the enclosing Coupling,
 *                used to resynchronize the whole coupling after any mesh
 *                mutation (see `Coupling::CollectAllVariables()`).
 */
template <class VAR, class PST>
void ProblemBase<VAR, PST>::save_amr(const std::vector<VAR*>& all_vars) {
  if (this->amr_ != nullptr) {
    this->amr_->StepDerefine(this->variables_, all_vars);
    this->amr_->StepRefine(this->variables_, all_vars);
  }
}

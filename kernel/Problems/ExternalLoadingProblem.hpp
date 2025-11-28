/**
 * @file ExternalLoadingProblem.hpp
 * @author cp273896 (clement.plumecocq@cea.fr)
 * @brief
 * @version 0.1
 * @date 2025-11-27
 *
 * Copyright CEA (C) 2025
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
 *
 */
#include <algorithm>
#include <cmath>
#include <functional>
#include <list>
#include <memory>
#include <mfem/linalg/vector.hpp>
#include <string>
#include <tuple>
#include <vector>

#include "Convergence/Convergence.hpp"
#include "Convergence/PhysicalConvergence.hpp"
#include "ExternalLoading/LoadingOptions.hpp"
#include "Parameters/Parameter.hpp"
#include "Problems/ProblemBase.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

template <template <CoordinateSystem> class LOADINGOPE, CoordinateSystem CS, class VAR, class PST>
class ExternalLoadingProblem {
 private:
  const Parameters& params_;
  std::vector<std::tuple<std::string, mfem::Vector>> get_coordinates();
  void get_parameters();

  std::shared_ptr<LOADINGOPE<CS>> LOPE;

 protected:
  std::string name_{"Unnamed problem"};

  VAR& variables_;
  std::vector<VAR*> auxvariables_;
  PST& pst_;
  std::vector<mfem::Vector> unknown_;

  std::vector<std::tuple<std::string, std::string>> sorted_chemical_system_;
  std::vector<std::tuple<std::string, bool, double>> var_convergence_;

 public:
  template <class... Args>
  ExternalLoadingProblem(const std::string& name, const Parameters& params, VAR& variables,
                         PST& pst, Args&&... auxvariable);

  /////////////////////////////////////////////////////
  void initialize(const double& initial_time);
  std::vector<std::tuple<std::string, bool, double>> get_convergence();
  /////////////////////////////////////////////////////
  void execute(const int& iter, double& next_time, const double& current_time,
               const double& current_time_step);
  void do_time_step(double& next_time, const double& current_time, double current_time_step,
                    const int iter, std::vector<mfem::Vector>& unks,
                    const std::vector<std::vector<std::string>>& unks_info);

  /////////////////////////////////////////////////////
  void save(const int& iter, const double& current_time);
  void post_execute(const int& iter, const double& current_time, const double& current_time_step);
  void post_processing(const int& iter, const double& current_time,
                       const double& current_time_step);
  void update();
  VAR get_problem_variables();
  std::string get_name();
  /////////////////////////////////////////////////////
  void finalize();
  /////////////////////////////////////////////////////

  ~ExternalLoadingProblem();
};

////////////////////////////////
////////////////////////////////

/**
 * @brief Construct a new External Loading Problem< L O A D I N G O P E,  C S,  V A R,  P S T>::
 * External Loading Problem object
 *
 * @tparam LOADINGOPE
 * @tparam CS
 * @tparam VAR
 * @tparam PST
 * @tparam Args
 * @param name
 * @param params
 * @param variables
 * @param pst
 * @param auxvariables
 */
template <template <CoordinateSystem> class LOADINGOPE, CoordinateSystem CS, class VAR, class PST>
template <class... Args>
ExternalLoadingProblem<LOADINGOPE, CS, VAR, PST>::ExternalLoadingProblem(const std::string& name,
                                                                         const Parameters& params,
                                                                         VAR& variables, PST& pst,
                                                                         Args&&... auxvariables)
    : name_(name), params_(params), variables_(variables), pst_(pst) {
  this->LOPE = std::make_shared<LOADINGOPE<CS>>(params_);
  if constexpr (sizeof...(auxvariables) == 0) {
    this->auxvariables_.resize(0);
  } else {
    this->auxvariables_ = {&auxvariables...};
  }
}

/**
 * @brief
 *
 * @tparam LOADINGOPE
 * @tparam CS
 * @tparam VAR
 * @tparam PST
 * @param initial_time
 */
template <template <CoordinateSystem> class LOADINGOPE, CoordinateSystem CS, class VAR, class PST>
void ExternalLoadingProblem<LOADINGOPE, CS, VAR, PST>::initialize(const double& initial_time) {
  this->LOPE->initialize();
}

/**
 * @brief
 *
 * @tparam LOADINGOPE
 * @tparam CS
 * @tparam VAR
 * @tparam PST
 */
template <template <CoordinateSystem> class LOADINGOPE, CoordinateSystem CS, class VAR, class PST>
void ExternalLoadingProblem<LOADINGOPE, CS, VAR, PST>::update() {
  for (std::size_t i = 0; i < this->variables_.get_variables_number(); i++) {
    auto& var = this->variables_.getIVariable(i);
    var.update(this->unknown_[i]);
  }
  this->unknown_.clear();
}

/**
 * @brief
 *
 * @tparam LOADINGOPE
 * @tparam CS
 * @tparam VAR
 * @tparam PST
 * @param iter
 * @param current_time
 * @param current_time_step
 */
template <template <CoordinateSystem> class LOADINGOPE, CoordinateSystem CS, class VAR, class PST>
void ExternalLoadingProblem<LOADINGOPE, CS, VAR, PST>::post_execute(
    const int& iter, const double& current_time, const double& current_time_step) {}

template <template <CoordinateSystem> class LOADINGOPE, CoordinateSystem CS, class VAR, class PST>
void ExternalLoadingProblem<LOADINGOPE, CS, VAR, PST>::do_time_step(
    double& next_time, const double& current_time, double current_time_step, const int iter,
    std::vector<mfem::Vector>& vect_unk, const std::vector<std::vector<std::string>>& unks_info) {
  next_time = current_time + current_time_step;
  std::vector<std::tuple<std::string, mfem::Vector>> coordinates_gf = this->get_coordinates();

  this->LOPE->do_time_step(next_time, this->unknown_, unks_info, coordinates_gf);
}

/**
 * @brief
 *
 * @tparam LOADINGOPE
 * @tparam CS
 * @tparam VAR
 * @tparam PST
 * @param iter
 * @param next_time
 * @param current_time
 * @param current_time_step
 */
template <template <CoordinateSystem> class LOADINGOPE, CoordinateSystem CS, class VAR, class PST>
void ExternalLoadingProblem<LOADINGOPE, CS, VAR, PST>::execute(const int& iter, double& next_time,
                                                               const double& current_time,
                                                               const double& current_time_step) {
  int rank = mfem::Mpi::WorldRank();
  if (rank == 0) {
    SlothInfo::verbose("   ============================== ");
    SlothInfo::verbose("   ==== Problem : ", this->name_);
    SlothInfo::verbose("   ============================== ");
  }

  std::vector<std::vector<std::string>> vect_unk_info;

  size_t num_vars = this->variables_.getVariables().size();
  unknown_.reserve(num_vars);
  vect_unk_info.reserve(num_vars);

  for (auto& var : this->variables_.getVariables()) {
    unknown_.push_back(var.get_unknown());

    auto& unk_info = vect_unk_info.emplace_back(var.get_additional_variable_info());
    unk_info.insert(unk_info.begin(), var.getVariableName());
  }

  this->do_time_step(next_time, current_time, current_time_step, iter, unknown_, vect_unk_info);
}

/**
 * @brief
 *
 * @tparam LOADINGOPE
 * @tparam CS
 * @tparam VAR
 * @tparam PST
 */
template <template <CoordinateSystem> class LOADINGOPE, CoordinateSystem CS, class VAR, class PST>
void ExternalLoadingProblem<LOADINGOPE, CS, VAR, PST>::get_parameters() {}

/**
 * @brief
 *
 * @tparam LOADINGOPE
 * @tparam CS
 * @tparam VAR
 * @tparam PST
 */
template <template <CoordinateSystem> class LOADINGOPE, CoordinateSystem CS, class VAR, class PST>
void ExternalLoadingProblem<LOADINGOPE, CS, VAR, PST>::finalize() {}

/**
 * @brief
 *
 * @tparam LOADINGOPE
 * @tparam CS
 * @tparam VAR
 * @tparam PST
 * @param iter
 * @param current_time
 * @param current_time_step
 */
template <template <CoordinateSystem> class LOADINGOPE, CoordinateSystem CS, class VAR, class PST>
void ExternalLoadingProblem<LOADINGOPE, CS, VAR, PST>::post_processing(
    const int& iter, const double& current_time, const double& current_time_step) {
  this->save(iter, current_time);
}

/**
 * @brief
 *
 * @tparam LOADINGOPE
 * @tparam CS
 * @tparam VAR
 * @tparam PST
 * @return std::vector<std::tuple<std::string, mfem::Vector>>
 */
template <template <CoordinateSystem> class LOADINGOPE, CoordinateSystem CS, class VAR, class PST>
std::vector<std::tuple<std::string, mfem::Vector>>
ExternalLoadingProblem<LOADINGOPE, CS, VAR, PST>::get_coordinates() {
  std::vector<std::tuple<std::string, mfem::Vector>> aux_gf;
  for (const auto& auxvar_vec : this->auxvariables_) {
    for (const auto& auxvar : auxvar_vec->getVariables()) {
      const auto gf = auxvar.get_unknown();
      auto variable_info = auxvar.get_additional_variable_info();
      const std::string& symbol = toUpperCase(auxvar.get_additional_variable_info().back());
      if ((symbol == "XCOORD") || (symbol == "YCOORD") || (symbol == "ZCOORD")) {
        MFEM_VERIFY(variable_info.size() == 1,
                    "Error while getting coordinates. Expected ['XYZCOORD']");
        aux_gf.emplace_back(std::make_tuple(symbol, gf));
      }
    }
  }

  return aux_gf;
}

/**
 * @brief
 *
 * @tparam LOADINGOPE
 * @tparam CS
 * @tparam VAR
 * @tparam PST
 * @return std::vector<std::tuple<std::string, bool, double>>
 */
template <template <CoordinateSystem> class LOADINGOPE, CoordinateSystem CS, class VAR, class PST>
std::vector<std::tuple<std::string, bool, double>>
ExternalLoadingProblem<LOADINGOPE, CS, VAR, PST>::get_convergence() {
  return this->var_convergence_;
}

/**
 * @brief
 *
 * @tparam LOADINGOPE
 * @tparam CS
 * @tparam VAR
 * @tparam PST
 * @return std::string
 */
template <template <CoordinateSystem> class LOADINGOPE, CoordinateSystem CS, class VAR, class PST>
std::string ExternalLoadingProblem<LOADINGOPE, CS, VAR, PST>::get_name() {
  return this->name_;
}

/**
 * @brief
 *
 * @tparam LOADINGOPE
 * @tparam CS
 * @tparam VAR
 * @tparam PST
 * @param iter
 * @param current_time
 */
template <template <CoordinateSystem> class LOADINGOPE, CoordinateSystem CS, class VAR, class PST>
void ExternalLoadingProblem<LOADINGOPE, CS, VAR, PST>::save(const int& iter,
                                                            const double& current_time) {
  auto vars = this->get_problem_variables();
  this->pst_.save_variables(vars, iter, current_time);
}

/**
 * @brief
 *
 * @tparam LOADINGOPE
 * @tparam CS
 * @tparam VAR
 * @tparam PST
 * @return VAR
 */
template <template <CoordinateSystem> class LOADINGOPE, CoordinateSystem CS, class VAR, class PST>
VAR ExternalLoadingProblem<LOADINGOPE, CS, VAR, PST>::get_problem_variables() {
  return this->variables_;
}

/**
 * @brief Destroy the External Loading Problem< L O A D I N G O P E,  C S,  V A R,  P S T>::
 * External Loading Problem object
 *
 * @tparam LOADINGOPE
 * @tparam CS
 * @tparam VAR
 * @tparam PST
 */
template <template <CoordinateSystem> class LOADINGOPE, CoordinateSystem CS, class VAR, class PST>
ExternalLoadingProblem<LOADINGOPE, CS, VAR, PST>::~ExternalLoadingProblem() {}

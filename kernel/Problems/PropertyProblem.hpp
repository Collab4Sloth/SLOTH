/**
 * @file PropertyProblem.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Class used to compute primary variables as a function of auxialiary variables
 * @version 0.1
 * @date 2025-05-09
 *
 * Copyright CEA (c) 2025
 *
 */
#include <algorithm>
#include <functional>
#include <list>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "Convergence/PhysicalConvergence.hpp"
#include "Parameters/Parameter.hpp"
#include "Problems/ProblemBase.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once
template <class PROPERTY, class VAR, class PST>
class Property_problem : public ProblemBase<VAR, PST> {
 private:
  PROPERTY* PP_;

 public:
  template <class... Args>
  Property_problem(const std::string& name, const Parameters& params, VAR& variables, PST& pst,
                   const PhysicalConvergence& convergence, Args&&... auxvariable);
  template <class... Args>
  Property_problem(const Parameters& params, VAR& variables, PST& pst,
                   const PhysicalConvergence& convergence, Args&&... auxvariable);

  void do_time_step(double& next_time, const double& current_time, double current_time_step,
                    const int iter, std::vector<std::unique_ptr<mfem::Vector>>& unks,
                    const std::vector<std::vector<std::string>>& unks_info) override;

  ~Property_problem();
};

////////////////////////////////
////////////////////////////////

/**
 * @brief Construct a new Property_problem<PROPERTY, VAR, PST>::Property_problem object
 *
 * @tparam PROPERTY
 * @tparam VAR
 * @tparam PST
 * @tparam Args
 * @param params
 * @param variables
 * @param pst
 * @param convergence
 * @param auxvariables
 */
template <class PROPERTY, class VAR, class PST>
template <class... Args>
Property_problem<PROPERTY, VAR, PST>::Property_problem(const Parameters& params, VAR& variables,
                                                       PST& pst,
                                                       const PhysicalConvergence& convergence,
                                                       Args&&... auxvariables)
    : ProblemBase<VAR, PST>("PropertyProblem", variables, pst, convergence, auxvariables...) {
  this->PP_ = new PROPERTY(params);
}

/**
 * @brief Construct a new Property_problem<PROPERTY, VAR, PST>::Property_problem object
 *
 * @tparam PROPERTY
 * @tparam VAR
 * @tparam PST
 * @tparam Args
 * @param name
 * @param params
 * @param variables
 * @param pst
 * @param convergence
 * @param auxvariables
 */
template <class PROPERTY, class VAR, class PST>
template <class... Args>
Property_problem<PROPERTY, VAR, PST>::Property_problem(const std::string& name,
                                                       const Parameters& params, VAR& variables,
                                                       PST& pst,
                                                       const PhysicalConvergence& convergence,
                                                       Args&&... auxvariables)
    : ProblemBase<VAR, PST>(name, variables, pst, convergence, auxvariables...) {
  this->PP_ = new PROPERTY(params);
}

/**
 * @brief  Do a time-step by calling Step method of the ODE

 *
 * @tparam PROPERTY
 * @tparam VAR
 * @tparam PST
 * @param next_time
 * @param current_time
 * @param current_time_step
 * @param iter
 * @param vect_unk
 * @param unks_info
 */
template <class PROPERTY, class VAR, class PST>
void Property_problem<PROPERTY, VAR, PST>::do_time_step(
    double& next_time, const double& current_time, double current_time_step, const int iter,
    std::vector<std::unique_ptr<mfem::Vector>>& vect_unk,
    const std::vector<std::vector<std::string>>& unks_info) {
  std::vector<mfem::Vector> tp_gf = this->get_tp_conditions();
  std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>
      output_system = this->get_output_system(unks_info, vect_unk);

  const size_t unk_size = vect_unk.size();

  std::vector<mfem::ParGridFunction> vect_aux_gf;
  std::vector<std::vector<std::string>> vect_aux_infos;
  for (const auto& auxvar_vec : this->auxvariables_) {
    for (const auto& auxvar : auxvar_vec->getVariables()) {
      // GF
      vect_aux_gf.emplace_back(std::move(auxvar.get_gf()));
      // Information
      std::vector<std::string> var_info = auxvar.get_additional_variable_info();
      var_info.push_back(auxvar.getVariableName());
      vect_aux_infos.emplace_back(std::move(var_info));
    }
  }

  this->PP_->compute(vect_unk, unks_info, vect_aux_gf, vect_aux_infos);

  // Recover unknowns
  for (size_t i = 0; i < unk_size; i++) {
    auto& unk_i = *(vect_unk[i]);
    this->unknown_.emplace_back(unk_i);
  }

  next_time = current_time + current_time_step;
}

/**
 * @brief Destroy the Property_problem<PROPERTY, VAR, PST>::Property_problem object
 *
 * @tparam PROPERTY
 * @tparam VAR
 * @tparam PST
 */
template <class PROPERTY, class VAR, class PST>
Property_problem<PROPERTY, VAR, PST>::~Property_problem() {}

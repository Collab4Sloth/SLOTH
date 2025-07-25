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

#include "Convergence/Convergence.hpp"
#include "Convergence/PhysicalConvergence.hpp"
#include "Parameters/Parameter.hpp"
#include "Problems/ProblemBase.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once
template <class PROPERTY, class VAR, class PST>
class Property_problem : public ProblemBase<VAR, PST> {
 private:
  PROPERTY* PP_;

  std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>
  get_output_system(const std::vector<std::vector<std::string>>& unks_info,
                    std::vector<std::unique_ptr<mfem::Vector>>& vect_unk);
  std::vector<std::tuple<std::vector<std::string>, mfem::Vector>> get_input_system();

 public:
  template <class... Args>
  Property_problem(const std::string& name, const Parameters& params, VAR& variables, PST& pst,
                   Convergence& convergence, Args&&... auxvariable);
  template <class... Args>
  Property_problem(const Parameters& params, VAR& variables, PST& pst, Convergence& convergence,
                   Args&&... auxvariable);

  template <class... Args>
  Property_problem(const std::string& name, const Parameters& params, VAR& variables, PST& pst,
                   Args&&... auxvariable);

  template <class... Args>
  Property_problem(const Parameters& params, VAR& variables, PST& pst, Args&&... auxvariable);

  void initialize(const double& initial_time) override;

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
                                                       PST& pst, Convergence& convergence,
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
 * @param params
 * @param variables
 * @param pst
 * @param auxvariables
 */
template <class PROPERTY, class VAR, class PST>
template <class... Args>
Property_problem<PROPERTY, VAR, PST>::Property_problem(const Parameters& params, VAR& variables,
                                                       PST& pst, Args&&... auxvariables)
    : ProblemBase<VAR, PST>("PropertyProblem", variables, pst, auxvariables...) {
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
                                                       PST& pst, Convergence& convergence,
                                                       Args&&... auxvariables)
    : ProblemBase<VAR, PST>(name, variables, pst, convergence, auxvariables...) {
  this->PP_ = new PROPERTY(params);
}

template <class PROPERTY, class VAR, class PST>
template <class... Args>
Property_problem<PROPERTY, VAR, PST>::Property_problem(const std::string& name,
                                                       const Parameters& params, VAR& variables,
                                                       PST& pst, Args&&... auxvariables)
    : ProblemBase<VAR, PST>(name, variables, pst, auxvariables...) {
  this->PP_ = new PROPERTY(params);
}

/**
 * @brief
 *
 * @tparam OPE
 * @tparam VAR
 * @tparam PST
 * @param initial_time
 */
template <class PROPERTY, class VAR, class PST>
void Property_problem<PROPERTY, VAR, PST>::initialize(const double& initial_time) {}

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
  // Primary variables
  std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>
      output_system = this->get_output_system(unks_info, vect_unk);
  // Auxiliary variables
  std::vector<std::tuple<std::vector<std::string>, mfem::Vector>> input_system =
      this->get_input_system();

  this->PP_->compute(output_system, input_system);

  // Recover unknowns
  const size_t unk_size = vect_unk.size();

  for (size_t i = 0; i < unk_size; i++) {
    auto& unk_i = *(vect_unk[i]);
    this->unknown_.emplace_back(unk_i);
  }

  next_time = current_time + current_time_step;
}

/**
 * @brief
 *
 * @tparam PROPERTY
 * @tparam VAR
 * @tparam PST
 * @param unks_info
 * @param vect_unk
 * @return std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>
 */
template <class PROPERTY, class VAR, class PST>
std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>
Property_problem<PROPERTY, VAR, PST>::get_output_system(
    const std::vector<std::vector<std::string>>& unks_info,
    std::vector<std::unique_ptr<mfem::Vector>>& vect_unk) {
  std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>
      output_system;
  for (size_t i = 0; i < vect_unk.size(); ++i) {
    output_system.emplace_back(unks_info[i], *vect_unk[i]);
  }

  return output_system;
}

/**
 * @brief
 *
 * @tparam PROPERTY
 * @tparam VAR
 * @tparam PST
 * @return std::vector<std::tuple<std::vector<std::string>, mfem::Vector>>
 */
template <class PROPERTY, class VAR, class PST>
std::vector<std::tuple<std::vector<std::string>, mfem::Vector>>
Property_problem<PROPERTY, VAR, PST>::get_input_system() {
  std::vector<std::tuple<std::vector<std::string>, mfem::Vector>> input_system;

  for (const auto& auxvar_vec : this->auxvariables_) {
    for (const auto& auxvar : auxvar_vec->getVariables()) {
      const auto gf = auxvar.get_unknown();

      auto variable_info = auxvar.get_additional_variable_info();
      input_system.emplace_back(std::make_tuple(variable_info, gf));
    }
  }

  return input_system;
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

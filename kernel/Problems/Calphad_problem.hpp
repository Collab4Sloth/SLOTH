/**
 * @file Calphad_problem.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Class used to defined a Calphad_Problem objet
 * @version 0.1
 * @date 2025-01-07
 *
 * Copyright CEA (c) 2025
 *
 */
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
template <class CALPHAD, class VAR, class PST>
class Calphad_Problem : public ProblemBase<VAR, PST> {
 private:
  CALPHAD* CC_;

  void check_variables_consistency();

  std::vector<mfem::Vector> get_tp_conditions();
  std::vector<std::tuple<std::string, std::string>> get_chemical_system();
  std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>
  get_output_system(const std::vector<std::vector<std::string>>& unks_info,
                    std::vector<std::unique_ptr<mfem::Vector>>& vect_unk);

 public:
  template <class... Args>
  Calphad_Problem(const Parameters& params, VAR& variables, PST& pst,
                  const PhysicalConvergence& convergence, std::list<int> pop_elem,
                  Args&&... auxvariable);

  template <class... Args>
  Calphad_Problem(const std::string& name, const Parameters& params, VAR& variables, PST& pst,
                  const PhysicalConvergence& convergence, std::list<int> pop_elem,
                  Args&&... auxvariable);

  template <class... Args>
  Calphad_Problem(const Parameters& params, VAR& variables, PST& pst,
                  const PhysicalConvergence& convergence, Args&&... auxvariable);

  template <class... Args>
  Calphad_Problem(const std::string& name, const Parameters& params, VAR& variables, PST& pst,
                  const PhysicalConvergence& convergence, Args&&... auxvariable);

  /////////////////////////////////////////////////////
  void initialize(const double& initial_time) override;

  /////////////////////////////////////////////////////

  void do_time_step(double& next_time, const double& current_time, double current_time_step,
                    const int iter, std::vector<std::unique_ptr<mfem::Vector>>& unks,
                    const std::vector<std::vector<std::string>>& unks_info) override;

  void get_parameters();

  /////////////////////////////////////////////////////
  void finalize() override;
  /////////////////////////////////////////////////////

  ~Calphad_Problem();
};

////////////////////////////////
////////////////////////////////

/**
 * @brief Construct a new  Calphad_Problem object (with auxiliary variables)
 * @warning At least two auxiliary variable. The first auxiliary variable is the temperature and the
 * seccond is the pressure.

 * @tparam CALPHAD
 * @tparam VAR
 * @tparam PST
 * @tparam Args
 * @param params
 * @param variables
 * @param pst
 * @param convergence
 * @param pop_elem
 * @param auxvariables
 */
template <class CALPHAD, class VAR, class PST>
template <class... Args>
Calphad_Problem<CALPHAD, VAR, PST>::Calphad_Problem(const Parameters& params, VAR& variables,
                                                    PST& pst,
                                                    const PhysicalConvergence& convergence,
                                                    std::list<int> pop_elem, Args&&... auxvariables)
    : ProblemBase<VAR, PST>("Calphad Problem", variables, pst, convergence, pop_elem,
                            auxvariables...) {
  this->check_variables_consistency();
  this->CC_ = new CALPHAD(params);
  this->get_parameters();
}

/**
 * @brief Construct a new Calphad_Problem<CALPHAD, VAR, PST>::Calphad_Problem object
 *
 * @tparam CALPHAD
 * @tparam VAR
 * @tparam PST
 * @tparam Args
 * @param name
 * @param params
 * @param variables
 * @param pst
 * @param convergence
 * @param pop_elem
 * @param auxvariables
 */
template <class CALPHAD, class VAR, class PST>
template <class... Args>
Calphad_Problem<CALPHAD, VAR, PST>::Calphad_Problem(const std::string& name,
                                                    const Parameters& params, VAR& variables,
                                                    PST& pst,
                                                    const PhysicalConvergence& convergence,
                                                    std::list<int> pop_elem, Args&&... auxvariables)
    : ProblemBase<VAR, PST>(name, variables, pst, convergence, pop_elem, auxvariables...) {
  this->check_variables_consistency();
  this->CC_ = new CALPHAD(params);
  this->get_parameters();
}

/**
 * @brief Construct a new Calphad_Problem<CALPHAD, VAR, PST>::Calphad_Problem object
 *
 * @tparam CALPHAD
 * @tparam VAR
 * @tparam PST
 * @tparam Args
 * @param params
 * @param variables
 * @param pst
 * @param convergence
 * @param auxvariables
 */
template <class CALPHAD, class VAR, class PST>
template <class... Args>
Calphad_Problem<CALPHAD, VAR, PST>::Calphad_Problem(const Parameters& params, VAR& variables,
                                                    PST& pst,
                                                    const PhysicalConvergence& convergence,
                                                    Args&&... auxvariables)
    : ProblemBase<VAR, PST>("Calphad Problem", variables, pst, convergence, auxvariables...) {
  this->check_variables_consistency();

  this->CC_ = new CALPHAD(params);
  this->get_parameters();
}

/**
 * @brief Construct a new Calphad_Problem<CALPHAD, VAR, PST>::Calphad_Problem object
 *
 * @tparam CALPHAD
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
template <class CALPHAD, class VAR, class PST>
template <class... Args>
Calphad_Problem<CALPHAD, VAR, PST>::Calphad_Problem(const std::string& name,
                                                    const Parameters& params, VAR& variables,
                                                    PST& pst,
                                                    const PhysicalConvergence& convergence,
                                                    Args&&... auxvariables)
    : ProblemBase<VAR, PST>(name, variables, pst, convergence, auxvariables...) {
  this->check_variables_consistency();
  this->CC_ = new CALPHAD(params);
  this->get_parameters();
}

/**
 * @brief Check consistency of variables. Name of a variable must have one of the following prefix:
 * xx_ for molar fraction,
 * nn_ for mole number,
 * mu_ for chemical potential,
 * df_ for driving force
 *
 */
template <class CALPHAD, class VAR, class PST>
void Calphad_Problem<CALPHAD, VAR, PST>::check_variables_consistency() {
  // ---------------------
  // Primary variables
  // ---------------------
  for (const auto& var : this->variables_.getVariables()) {
    MFEM_VERIFY(var.get_additional_variable_info().size() > 0,
                "Calphad problems requires ouputs with at least two additional informations  \n");

    switch (calphad_outputs::from(var.get_additional_variable_info().back())) {
      case calphad_outputs::mu: {
        MFEM_VERIFY(var.get_additional_variable_info().size() == 2,
                    "Calphad problems requires that chemical potential ouputs are defined with two "
                    "additional informations: first an element, second the symbol 'mu'. \n");

        SlothInfo::debug("Output : chemical potential for ", var.get_additional_variable_info()[0]);
        break;
      }
      case calphad_outputs::x: {
        MFEM_VERIFY(var.get_additional_variable_info().size() == 3,
                    "Calphad problems requires that element molar fraction ouputs are defined with "
                    "three additional informations: first an element, second a phase and third "
                    "the symbol 'x'. \n");
        SlothInfo::debug("Output : molar fraction for ", var.get_additional_variable_info()[0]);
        break;
      }
      case calphad_outputs::g:
      case calphad_outputs::gm: {
        MFEM_VERIFY(var.get_additional_variable_info().size() == 2,
                    "Calphad problems requires that Gibbs energy ouputs are defined with two "
                    "additional informations: first a phase,  second the symbol 'g' or 'gm'. \n");
        SlothInfo::debug("Output : Gibbs energy for ", var.get_additional_variable_info()[0]);
        break;
      }
      case calphad_outputs::h:
      case calphad_outputs::hm: {
        MFEM_VERIFY(var.get_additional_variable_info().size() == 2,
                    "Calphad problems requires that enthalpy ouputs are defined with two "
                    "additional informations: first a phase,  second the symbol  'h' or 'hm'. \n");
        SlothInfo::debug("Output : enthalpy for ", var.get_additional_variable_info()[0]);
        break;
      }
      case calphad_outputs::dgm: {
        MFEM_VERIFY(var.get_additional_variable_info().size() == 2,
                    "Calphad problems requires that driving forces ouputs are defined with the "
                    "name of the phase and the symbol 'dgm' as the additional information. \n");
        SlothInfo::debug("Output : driving force for ", var.get_additional_variable_info()[0]);
        break;
      }
      case calphad_outputs::mob: {
        MFEM_VERIFY(
            var.get_additional_variable_info().size() == 3,
            "Calphad problems requires that mobility ouputs are defined with two "
            "additional informations: first a phase, second an element, the symbol 'mob'. \n");

        SlothInfo::debug("Output : mobility for ", var.get_additional_variable_info()[0]);
        break;
      }
      case calphad_outputs::cp: {
        MFEM_VERIFY(var.get_additional_variable_info().size() == 1,
                    "Calphad problems requires that heat capacity ouputs are defined with the "
                    "symbol cp as the only additional information. \n");
        SlothInfo::debug("Output : heat capacity ");
        break;
      }
    }
  }
  // ---------------------
  // Auxiliary variables
  // ---------------------
  MFEM_VERIFY(this->auxvariables_.size() > 2,
              "Calphad problems requires three set of auxiliary variables: the temperature, the "
              "pressure and the composition");
  // Temperature and Pressure
  bool has_temperature = false;
  bool has_pressure = false;
  for (const auto& auxvar_vec : this->auxvariables_) {
    for (const auto& auxvar : auxvar_vec->getVariables()) {
      if (auxvar.get_additional_variable_info()[0] == "Temperature") {
        has_temperature = true;
      }
      if (auxvar.get_additional_variable_info()[0] == "Pressure") {
        has_pressure = true;
      }
    }
  }
  MFEM_VERIFY(has_temperature && has_pressure,
              "Calphad problems requires  temperature and pressure as auxiliary variables");

  // Composition
  // TODO(cci) : add check of the type of variable among allowed type
}

/**
 * @brief Initialization of the problem depending on the initialization of the CalphadBase objet
 *
 * @tparam CALPHAD
 * @tparam VAR
 * @tparam PST
 * @param initial_time
 */
template <class CALPHAD, class VAR, class PST>
void Calphad_Problem<CALPHAD, VAR, PST>::initialize(const double& initial_time) {
  this->CC_->initialize();
}

/**
 * @brief  Do a time-step by calling Step method of the ODE

 *
 * @tparam CALPHAD
 * @tparam VAR
 * @tparam PST
 * @param next_time
 * @param current_time
 * @param current_time_step
 * @param iter
 * @param vect_unk
 * @param unks_info
 */
template <class CALPHAD, class VAR, class PST>
void Calphad_Problem<CALPHAD, VAR, PST>::do_time_step(
    double& next_time, const double& current_time, double current_time_step, const int iter,
    std::vector<std::unique_ptr<mfem::Vector>>& vect_unk,
    const std::vector<std::vector<std::string>>& unks_info) {
  int rank = mfem::Mpi::WorldRank();

  // Temperature, Pressure and composition comes from auxiliary variables
  // They can be updated at each time-step by other problem (thermal, inter-diffusion,...) or stay
  // constant

  std::vector<mfem::Vector> tp_gf = this->get_tp_conditions();
  std::vector<std::tuple<std::string, std::string>> chemical_system = this->get_chemical_system();
  std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>
      output_system = this->get_output_system(unks_info, vect_unk);

  const size_t unk_size = vect_unk.size();

  this->CC_->execute(iter, tp_gf, chemical_system, output_system);

  for (size_t i = 0; i < unk_size; i++) {
    auto& unk_i = *(vect_unk[i]);
    this->unknown_.emplace_back(unk_i);
  }

  next_time = current_time + current_time_step;
}

/**
 * @brief Get all parameters associated with the problem.
 *
 * @tparam CALPHAD
 * @tparam VAR
 * @tparam PST
 */
template <class CALPHAD, class VAR, class PST>
void Calphad_Problem<CALPHAD, VAR, PST>::get_parameters() {
  this->CC_->get_parameters();
}

/**
 * @brief Finalization of the problem depending on the finalization of the CalphadBase objet
 *
 * @tparam CALPHAD
 * @tparam VAR
 * @tparam PST
 */
template <class CALPHAD, class VAR, class PST>
void Calphad_Problem<CALPHAD, VAR, PST>::finalize() {
  this->CC_->finalize();
}

/**
 * @brief Define the temperature and the pressure use to calculate  equilibria
 * @warning By convention temperature, pressure are given in the first Variables objet Composition
 * is  given in the second Variables objet
 * @tparam CALPHAD
 * @tparam VAR
 * @tparam PST
 * @return std::vector<mfem::Vector>
 */
template <class CALPHAD, class VAR, class PST>
std::vector<mfem::Vector> Calphad_Problem<CALPHAD, VAR, PST>::get_tp_conditions() {
  std::vector<mfem::Vector> aux_gf;

  for (const auto& auxvar_vec : this->auxvariables_) {
    for (const auto& auxvar : auxvar_vec->getVariables()) {
      const auto gf = auxvar.get_unknown();
      aux_gf.emplace_back(gf);
    }
  }

  return aux_gf;
}

/**
 * @brief Define the composition use to calculate  equilibria
 * @warning By convention temperature, pressure are given in the first Variables objet Composition
 * is  given in the second Variables objet
 *
 * @tparam CALPHAD
 * @tparam VAR
 * @tparam PST
 * @return std::vector<std::tuple<std::string, std::string>>
 */
template <class CALPHAD, class VAR, class PST>
std::vector<std::tuple<std::string, std::string>>
Calphad_Problem<CALPHAD, VAR, PST>::get_chemical_system() {
  std::vector<std::tuple<std::string, std::string>> chemical_system;

  for (const auto& auxvar_vec : this->auxvariables_) {
    for (const auto& auxvar : auxvar_vec->getVariables()) {
      auto variable_info = auxvar.get_additional_variable_info();
      // Check consistency of additional info
      // For composition, 2 information are required : element symbol, unit symbol
      if (variable_info[0] == "Temperature") continue;
      if (variable_info[0] == "Pressure") continue;
      MFEM_VERIFY(
          variable_info.size() == 2,
          "Error while getting chemical system. Three  additional informations are excepted "
          ": the element symbol, the unit symbol");
      chemical_system.emplace_back(std::make_tuple(variable_info[0], variable_info[1]));
    }
  }
  return chemical_system;
}

/**
 * @brief Define the output system
 *
 * @tparam CALPHAD
 * @tparam VAR
 * @tparam PST
 * @param unks_info
 * @param vect_unk
 * @return std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>
 */
template <class CALPHAD, class VAR, class PST>
std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>
Calphad_Problem<CALPHAD, VAR, PST>::get_output_system(
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
 * @brief Destroy the Calphad_Problem<CALPHAD, VAR, PST>::Calphad_Problem object
 *
 * @tparam CALPHAD
 * @tparam VAR
 * @tparam PST
 */
template <class CALPHAD, class VAR, class PST>
Calphad_Problem<CALPHAD, VAR, PST>::~Calphad_Problem() {}

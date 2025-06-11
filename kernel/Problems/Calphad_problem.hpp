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
template <class CALPHAD, class VAR, class PST>
class Calphad_Problem : public ProblemBase<VAR, PST> {
 private:
  const Parameters& params_;
  bool is_KKS_;
  CALPHAD* CC_;
  int findIndexOfTuple(const std::vector<std::tuple<std::string, std::string>>& vec,
                       const std::string& target);
  void check_variables_consistency();

  std::vector<mfem::Vector> get_tp_conditions();
  std::vector<mfem::Vector> get_old_tp_conditions();
  std::tuple<std::string, mfem::Vector, mfem::Vector> get_phasefields();

  std::vector<std::tuple<std::string, std::vector<double>>> get_coord();

  std::vector<std::tuple<std::string, std::string, mfem::Vector, mfem::Vector>>
  get_molar_fractions();
  std::vector<std::tuple<std::string, std::string>> get_chemical_system();
  void check_phasefield();
  void check_molar_fractions();

  std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>
  get_output_system(const std::vector<std::vector<std::string>>& unks_info,
                    std::vector<std::unique_ptr<mfem::Vector>>& vect_unk);
  std::vector<std::tuple<std::vector<std::string>, mfem::Vector>> get_previous_output_system();

  void get_parameters();

 protected:
  std::vector<std::tuple<std::string, std::string>> sorted_chemical_system_;
  std::vector<std::string> sorted_KKS_phases_;

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
                            auxvariables...),
      params_(params) {
  // Mandatory to be placed before CALPHAD pointer creation
  this->get_parameters();
  this->check_variables_consistency();

  this->CC_ = new CALPHAD(params, this->is_KKS_);

  this->sorted_chemical_system_ = this->get_chemical_system();
  if (this->is_KKS_) {
    this->check_phasefield();
    this->check_molar_fractions();
  }
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
    : ProblemBase<VAR, PST>(name, variables, pst, convergence, pop_elem, auxvariables...),
      params_(params) {
  // Mandatory to be placed before CALPHAD pointer creation
  this->get_parameters();
  this->check_variables_consistency();

  this->CC_ = new CALPHAD(params, this->is_KKS_);
  this->sorted_chemical_system_ = this->get_chemical_system();
  if (this->is_KKS_) {
    this->check_phasefield();
    this->check_molar_fractions();
  }
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
    : ProblemBase<VAR, PST>("Calphad Problem", variables, pst, convergence, auxvariables...),
      params_(params) {
  // Mandatory to be placed before CALPHAD pointer creation
  this->get_parameters();
  this->check_variables_consistency();

  this->CC_ = new CALPHAD(params, this->is_KKS_);
  this->sorted_chemical_system_ = this->get_chemical_system();
  if (this->is_KKS_) {
    this->check_phasefield();
    this->check_molar_fractions();
  }
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
    : ProblemBase<VAR, PST>(name, variables, pst, convergence, auxvariables...), params_(params) {
  // Mandatory to be placed before CALPHAD pointer creation
  this->get_parameters();
  this->check_variables_consistency();
  this->CC_ = new CALPHAD(params, this->is_KKS_);
  this->sorted_chemical_system_ = this->get_chemical_system();
  if (this->is_KKS_) {
    this->check_phasefield();
    this->check_molar_fractions();
  }
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

    const std::string& symbol = toLowerCase(var.get_additional_variable_info().back());
    switch (calphad_outputs::from(symbol)) {
      case calphad_outputs::mu: {
        MFEM_VERIFY(var.get_additional_variable_info().size() == 2,
                    "Calphad problems requires that chemical potential ouputs are defined with two "
                    "additional informations: first an element, second the symbol 'mu'. \n");

        SlothInfo::debug("Output : chemical potential for ", var.get_additional_variable_info()[0]);
        break;
      }
      case calphad_outputs::dmu: {
        MFEM_VERIFY(var.get_additional_variable_info().size() == 2,
                    "Calphad problems requires that diffusion chemical potential ouputs are "
                    "defined with two "
                    "additional informations: first an element and second the symbol 'dmu'. \n");

        SlothInfo::debug("Output : diffusion chemical potential for ",
                         var.get_additional_variable_info()[0]);
        break;
      }
      case calphad_outputs::xph: {
        MFEM_VERIFY(var.get_additional_variable_info().size() == 2,
                    "Calphad problems requires that phase molar fraction ouputs are defined with "
                    "two additional informations: first a phase and second the symbol 'xph'. \n");
        SlothInfo::debug("Output : molar fraction for ", var.get_additional_variable_info()[0]);
        break;
      }
      case calphad_outputs::xp: {
        MFEM_VERIFY(var.get_additional_variable_info().size() == 3,
                    "Calphad problems requires that element molar fraction ouputs are defined with "
                    "three additional informations: first an element, second a phase and third "
                    "the symbol 'xp'. \n");
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
      const std::string& symbol = toUpperCase(auxvar.get_additional_variable_info().back());
      if (symbol == "T") {
        has_temperature = true;
      }
      if (symbol == "P") {
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
  this->CC_->initialize(this->sorted_chemical_system_);
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
  std::vector<mfem::Vector> tp_gf = this->get_tp_conditions();
  std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>
      output_system = this->get_output_system(unks_info, vect_unk);

  std::vector<std::tuple<std::vector<std::string>, mfem::Vector>> previous_output =
      get_previous_output_system();

  const size_t unk_size = vect_unk.size();

  // Execute
  if (!this->is_KKS_) {
    this->CC_->global_execute(iter, current_time_step, tp_gf, this->sorted_chemical_system_,
                              output_system, previous_output);
  } else {
    // Specific treatment for KKS problems
    std::tuple<std::string, mfem::Vector, mfem::Vector> phasefields_gf = this->get_phasefields();
    std::vector<std::tuple<std::string, std::vector<double>>> coordinates = this->get_coord();

    std::vector<mfem::Vector> tp_gf_old = this->get_old_tp_conditions();
    // vector<element, phase, x, x_old>
    std::vector<std::tuple<std::string, std::string, mfem::Vector, mfem::Vector>> x_phase_gf =
        this->get_molar_fractions();

    this->CC_->global_execute(iter, current_time_step, tp_gf, this->sorted_chemical_system_,
                              output_system, previous_output, phasefields_gf, tp_gf_old, x_phase_gf,
                              coordinates);
  }

  // Recover unknowns
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
  this->is_KKS_ = this->params_.template get_param_value_or_default<bool>("enable_KKS", false);
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
 * @brief Return unknown associated with Phase-fields at current and previous time-step
 *
 * @tparam CALPHAD
 * @tparam VAR
 * @tparam PST
 * @return std::vector<mfem::Vector>
 */
template <class CALPHAD, class VAR, class PST>
std::tuple<std::string, mfem::Vector, mfem::Vector>
Calphad_Problem<CALPHAD, VAR, PST>::get_phasefields() {
  std::tuple<std::string, mfem::Vector, mfem::Vector> aux_gf;

  for (const auto& auxvar_vec : this->auxvariables_) {
    for (const auto& auxvar : auxvar_vec->getVariables()) {
      const auto gf = auxvar.get_unknown();
      const auto gf_old = auxvar.get_second_to_last();
      auto variable_info = auxvar.get_additional_variable_info();
      const std::string& symbol = toLowerCase(auxvar.get_additional_variable_info().back());

      if (symbol != "phi") continue;
      MFEM_VERIFY(variable_info.size() == 2,
                  "Error while getting phase_field. Expected [name of the phase, 'phi']");
      aux_gf = std::make_tuple(variable_info[0], gf, gf_old);
    }
  }

  return aux_gf;
}

/**
 * @brief Return the coordinates as {(X1,{x1,....,xn}),(X2,{y1,....,yn}),(X3,{z1,....,zn})}
 *
 * @tparam CALPHAD
 * @tparam VAR
 * @tparam PST
 * @return std::vector<std::tuple<std::string, std::vector<double>>>
 */
template <class CALPHAD, class VAR, class PST>
std::vector<std::tuple<std::string, std::vector<double>>>
Calphad_Problem<CALPHAD, VAR, PST>::get_coord() {
  return this->variables_.getIVariable(0).get_coordinates();
}

/**
 * @brief Return unknown associated with element molar fractions by phase at current and previous
 * time-step (vector<element, phase, x, x_old>)
 * @remark Element are sorted as the chemical system read for equilibrium calculations
 *
 * @tparam CALPHAD
 * @tparam VAR
 * @tparam PST
 * @return std::vector<std::tuple<std::string,std::string, mfem::Vector, mfem::Vector>>
 */
template <class CALPHAD, class VAR, class PST>
std::vector<std::tuple<std::string, std::string, mfem::Vector, mfem::Vector>>
Calphad_Problem<CALPHAD, VAR, PST>::get_molar_fractions() {
  const auto size_v = this->sorted_chemical_system_.size();
  std::vector<std::tuple<std::string, std::string, mfem::Vector, mfem::Vector>> aux_gf;
  for (const auto& var : this->variables_.getVariables()) {
    const std::string& symbol = toLowerCase(var.get_additional_variable_info().back());
    if (var.get_additional_variable_info().size() == 3 &&
        (calphad_outputs::from(symbol) == calphad_outputs::xp)) {
      const auto gf = var.get_unknown();
      const auto gf_old = var.get_second_to_last();
      auto variable_info = var.get_additional_variable_info();
      // const int id = this->findIndexOfTuple(this->sorted_chemical_system_, variable_info[0]);

      aux_gf.emplace_back(std::make_tuple(variable_info[0], variable_info[1], gf, gf_old));
    }
  }

  return aux_gf;
}

template <class CALPHAD, class VAR, class PST>
std::vector<std::tuple<std::vector<std::string>, mfem::Vector>>
Calphad_Problem<CALPHAD, VAR, PST>::get_previous_output_system() {
  std::vector<std::tuple<std::vector<std::string>, mfem::Vector>> aux_gf;
  for (const auto& var : this->variables_.getVariables()) {
    const auto gf = var.get_unknown();

    aux_gf.emplace_back(std::make_tuple(var.get_additional_variable_info(), gf));
  }

  return aux_gf;
}

/**
 * @brief Get temperature, pressure and composition used to calculate  equilibria
 * @warning By convention temperature, pressure are given in the first Variables objet, Composition
 * is given in the second Variables objet
 * @tparam CALPHAD
 * @tparam VAR
 * @tparam PST
 * @return std::vector<mfem::Vector>
 */
template <class CALPHAD, class VAR, class PST>
std::vector<mfem::Vector> Calphad_Problem<CALPHAD, VAR, PST>::get_tp_conditions() {
  const auto size_v = this->sorted_chemical_system_.size() + 2;
  std::vector<mfem::Vector> aux_gf(size_v);

  for (const auto& auxvar_vec : this->auxvariables_) {
    for (const auto& auxvar : auxvar_vec->getVariables()) {
      const auto gf = auxvar.get_unknown();
      auto variable_info = auxvar.get_additional_variable_info();
      MFEM_VERIFY(!variable_info.empty(), "Empty variable_info encountered.");
      const std::string& symbol = toUpperCase(variable_info.back());
      if (symbol == "PHI") continue;
      if (symbol == "XCOORD") continue;
      if (symbol == "YCOORD") continue;
      if (symbol == "T") {
        aux_gf[0] = gf;
      } else if (symbol == "P") {
        aux_gf[1] = gf;
      } else if (symbol == "X" || symbol == "N") {
        MFEM_VERIFY(variable_info.size() == 2,
                    " Calphad_Problem<CALPHAD, VAR, PST>::get_tp_conditions() : expected "
                    "[component, 'x'] or [component, 'N']");
        const int id = this->findIndexOfTuple(this->sorted_chemical_system_, variable_info[0]);

        aux_gf[id + 2] = gf;
      }
    }
  }

  return aux_gf;
}

/**
 * @brief Get the values at point expansion for temperature, pressure and composition used to
 * calculate equilibria in KKS model
 * @warning By convention temperature, pressure are given in the first Variables objet, Composition
 * is given in the second Variables objet
 * @tparam CALPHAD
 * @tparam VAR
 * @tparam PST
 * @return std::vector<mfem::Vector>
 */
template <class CALPHAD, class VAR, class PST>
std::vector<mfem::Vector> Calphad_Problem<CALPHAD, VAR, PST>::get_old_tp_conditions() {
  const auto size_v = this->sorted_chemical_system_.size() + 2;
  std::vector<mfem::Vector> aux_gf(size_v);

  for (const auto& auxvar_vec : this->auxvariables_) {
    for (const auto& auxvar : auxvar_vec->getVariables()) {
      // Get values at previous time-step (see Variable.hpp)
      const auto gf = auxvar.get_second_to_last();
      //
      auto variable_info = auxvar.get_additional_variable_info();
      MFEM_VERIFY(!variable_info.empty(), "Empty variable_info encountered.");
      const std::string& symbol = toUpperCase(variable_info.back());
      if (symbol == "PHI") continue;
      if (symbol == "XCOORD") continue;
      if (symbol == "YCOORD") continue;
      if (symbol == "T") {
        aux_gf[0] = gf;
      } else if (symbol == "P") {
        aux_gf[1] = gf;
      } else if (symbol == "X" || symbol == "N") {
        MFEM_VERIFY(variable_info.size() == 2,
                    " Calphad_Problem<CALPHAD, VAR, PST>::get_tp_conditions() : expected "
                    "[component, 'x'] or [component, 'N']");
        const int id =
            this->findIndexOfTuple(this->sorted_chemical_system_, toUpperCase(variable_info[0]));

        aux_gf[id + 2] = gf;
      }
    }
  }

  return aux_gf;
}

/**
 * @brief Find index of first element of a tuple inside a vector of tuple
 *
 * @tparam CALPHAD
 * @tparam VAR
 * @tparam PST
 * @param vec
 * @param target
 * @return int
 */
template <class CALPHAD, class VAR, class PST>
int Calphad_Problem<CALPHAD, VAR, PST>::findIndexOfTuple(
    const std::vector<std::tuple<std::string, std::string>>& vec, const std::string& target) {
  auto it = std::find_if(vec.begin(), vec.end(),
                         [&target](const std::tuple<std::string, std::string>& t) {
                           return std::get<0>(t) == target;
                         });

  if (it != vec.end()) {
    return std::distance(vec.begin(), it);
  } else {
    return -1;  // Return -1 if the element is not found
  }
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
      MFEM_VERIFY(!variable_info.empty(), "Empty variable_info encountered.");
      const std::string& symbol = toUpperCase(variable_info.back());
      // Check consistency of additional info
      // For composition, 2 information are required : element symbol, unit symbol
      if (symbol == "T") continue;
      if (symbol == "P") continue;
      if (symbol == "PHI") continue;
      if (symbol == "XCOORD") continue;
      if (symbol == "YCOORD") continue;
      MFEM_VERIFY(variable_info.size() == 2,
                  "Error while getting chemical system. Two additional informations are excepted "
                  ": the element symbol, the unit symbol");
      chemical_system.emplace_back(std::make_tuple(toUpperCase(variable_info[0]), symbol));
    }
  }
  // Sort by alphabetical order
  std::ranges::sort(chemical_system,
                    [](const auto& a, const auto& b) { return std::get<0>(a) < std::get<0>(b); });

  return chemical_system;
}

/**
 * @brief Return the list of phases involved in the KKS problem
 *
 * @tparam CALPHAD
 * @tparam VAR
 * @tparam PST
 * @return std::vector<std::string>
 */
template <class CALPHAD, class VAR, class PST>
void Calphad_Problem<CALPHAD, VAR, PST>::check_phasefield() {
  bool has_phasefield = false;

  for (const auto& auxvar_vec : this->auxvariables_) {
    for (const auto& auxvar : auxvar_vec->getVariables()) {
      auto variable_info = auxvar.get_additional_variable_info();
      const std::string& symbol = toUpperCase(variable_info.back());
      // Check consistency of additional info
      // For phase-fields, 2 information are required : "PhaseField", the name of the phase
      if (symbol != "PHI") continue;
      MFEM_VERIFY(
          variable_info.size() == 2,
          "Error while getting phases for KKS problems. Two additional informations are excepted "
          ": PhaseField and the name of the phase");

      has_phasefield = true;
    }
  }
  MFEM_VERIFY(has_phasefield,
              "Error while getting phases for KKS problems. One phase-field variable is expected "
              "as auxiliary variable");
}

template <class CALPHAD, class VAR, class PST>
void Calphad_Problem<CALPHAD, VAR, PST>::check_molar_fractions() {
  // Molar fractions by phase for KKS problem
  // TODO(cci ): check if well defined
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

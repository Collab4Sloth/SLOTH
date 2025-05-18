/**
 * @file SlothNLFormIntegrator.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Class prodiving common methods to all NLFormIntegrators
 * @version 0.1
 * @date 2025-01-22
 *
 * @copyright Copyright (c) 2025
 *
 */
#include <algorithm>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Profiling/Profiling.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

/**
 * @brief Class prodiving common methods to all NLFormIntegrators
 *
 */
template <class VARS>
class SlothNLFormIntegrator {
 private:
  void manage_auxiliary_variables(std::vector<VARS*> auxvars);
  std::vector<mfem::ParGridFunction> vect_aux_gf_;
  std::vector<mfem::Vector> vect_aux_old_gf_;
  std::vector<std::vector<std::string>> vect_aux_infos_;

 protected:
  std::vector<VARS*> auxvariables_;
  const Parameters params_;

  std::vector<mfem::ParGridFunction> get_aux_gf();
  std::vector<mfem::Vector> get_aux_old_gf();
  std::vector<std::vector<std::string>> get_aux_infos();

 public:
  SlothNLFormIntegrator(const Parameters& params, std::vector<VARS*> auxvars);
  ~SlothNLFormIntegrator();
};
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Construct a new SlothNLFormIntegrator::SlothNLFormIntegrator object
 *
 * @tparam VARS
 * @param params
 * @param auxvars
 */
template <class VARS>
SlothNLFormIntegrator<VARS>::SlothNLFormIntegrator(const Parameters& params,
                                                   std::vector<VARS*> auxvars)
    : params_(params) {
  this->manage_auxiliary_variables(auxvars);
}

/**
 * @brief Return the vector of grid functions associated with the auxiliary variables
 *
 * @tparam VARS
 * @return std::vector<mfem::ParGridFunction>
 */
template <class VARS>
void SlothNLFormIntegrator<VARS>::manage_auxiliary_variables(std::vector<VARS*> auxvars) {
  this->auxvariables_ = auxvars;
  for (const auto& auxvar_vec : this->auxvariables_) {
    for (const auto& auxvar : auxvar_vec->getVariables()) {
      // GF
      this->vect_aux_gf_.emplace_back(std::move(auxvar.get_gf()));
      // GF at previous time-step
      this->vect_aux_old_gf_.emplace_back(std::move(auxvar.get_second_to_last()));

      // Information
      std::vector<std::string> var_info = auxvar.get_additional_variable_info();
      var_info.push_back(auxvar.getVariableName());

      this->vect_aux_infos_.emplace_back(std::move(var_info));
    }
  }
}

/**
 * @brief Return a vector of GridFunction associated with auxiliary variables
 * @remark Order of the vector is implicitly the same as the order of auxiliary variables
 *
 * @tparam VARS
 * @return std::vector<mfem::ParGridFunction>
 */
template <class VARS>
std::vector<mfem::ParGridFunction> SlothNLFormIntegrator<VARS>::get_aux_gf() {
  return vect_aux_gf_;
}

/**
 * @brief Return a vector of GridFunction associated with auxiliary variables at the previous
 * time-step
 * @remark Order of the vector is implicitly the same as the order of auxiliary variables
 *
 * @tparam VARS
 * @return std::vector<mfem::ParGridFunction>
 */
template <class VARS>
std::vector<mfem::Vector> SlothNLFormIntegrator<VARS>::get_aux_old_gf() {
  return vect_aux_old_gf_;
}

/**
 * @brief Return a vector of the additional information (vector of string) associated with auxiliary
 * variables.
 * @remark Order of the vector is implicitly the same as the order of auxiliary variables
 * @remark The last element of the vector is the name of the variable
 *
 * @tparam VARS
 * @return std::vector<std::vector<std::string>>
 */
template <class VARS>
std::vector<std::vector<std::string>> SlothNLFormIntegrator<VARS>::get_aux_infos() {
  return this->vect_aux_infos_;
}

/**
 * @brief Destroy the SlothNLFormIntegrator::SlothNLFormIntegrator object
 *
 * @tparam VARS
 */
template <class VARS>
SlothNLFormIntegrator<VARS>::~SlothNLFormIntegrator() {}

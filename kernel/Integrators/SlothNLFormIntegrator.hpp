/**
 * @file SlothNLFormIntegrator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Class prodiving common methods to all NLFormIntegrators
 * @version 0.1
 * @date 2025-09-05
 *
 * @anchor SlothNLFormIntegrator
 *
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
 */
#include <algorithm>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "Coefficients/Coefficient.hpp"
#include "Coefficients/Coefficients.hpp"
#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

/**
 * @brief Class prodiving common methods to all NLFormIntegrators
 *
 */
template <class VARS>
class SlothNLFormIntegrator : public mfem::BlockNonlinearFormIntegrator {
 private:
  void manage_auxiliary_variables(std::vector<VARS*> auxvars);
  std::vector<mfem::ParGridFunction> vect_aux_gf_;
  std::vector<mfem::Vector> vect_aux_old_gf_;
  std::vector<std::vector<std::string>> vect_aux_infos_;

 protected:
  std::string integrator_name_ = "";
  virtual void AssembleElementVector(const mfem::Array<const mfem::FiniteElement*>& el,
                                     mfem::ElementTransformation& Tr,
                                     const mfem::Array<const mfem::Vector*>& elfun,
                                     const mfem::Array<mfem::Vector*>& elvect) = 0;

  virtual void AssembleElementGrad(const mfem::Array<const mfem::FiniteElement*>& el,
                                   mfem::ElementTransformation& Tr,
                                   const mfem::Array<const mfem::Vector*>& elfun,
                                   const mfem::Array2D<mfem::DenseMatrix*>& elmat) = 0;
  virtual void get_coefficients() = 0;

  std::vector<mfem::ParGridFunction> u_old_;
  std::vector<mfem::ParGridFunction> aux_gf_;
  std::vector<mfem::Vector> aux_old_gf_;
  std::vector<std::vector<std::string>> aux_infos_;

  std::vector<VARS*> auxvariables_;
  Parameters params_;
  std::vector<Coefficients> coefficients_;
  unsigned int nb_blk_;

  std::vector<mfem::ParGridFunction> get_aux_gf();
  std::vector<mfem::Vector> get_aux_old_gf();
  std::vector<std::vector<std::string>> get_aux_infos();

  void check_coefficient_types(std::list<GlossaryType> expected_type);
  std::optional<Coefficient> get_coefficient(const int blk, GlossaryType type, unsigned int id);

 public:
  virtual void init() = 0;
  SlothNLFormIntegrator(const std::vector<mfem::ParGridFunction> u_old, const Parameters& params,
                        std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients);
  virtual ~SlothNLFormIntegrator() = default;
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
SlothNLFormIntegrator<VARS>::SlothNLFormIntegrator(const std::vector<mfem::ParGridFunction> u_old,
                                                   const Parameters& params,
                                                   std::vector<VARS*> auxvars,
                                                   const std::vector<Coefficients>& coefficients)
    : u_old_(u_old), params_(params), coefficients_(coefficients) {
  this->nb_blk_ = this->u_old_.size();
  this->manage_auxiliary_variables(auxvars);
  this->aux_gf_ = this->get_aux_gf();
  this->aux_old_gf_ = this->get_aux_old_gf();
  this->aux_infos_ = this->get_aux_infos();
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
      // var_info.push_back(auxvar.getVariableName());

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
 * @brief Check if the Coefficients passed to the SLOTH-based Integrator are consistent with its
 * expected list of GlossaryType
 *
 * @tparam VARS
 * @param expected_types
 */
template <class VARS>
void SlothNLFormIntegrator<VARS>::check_coefficient_types(std::list<GlossaryType> expected_types) {
  for (auto coefficients : this->coefficients_) {
    auto vect_types = coefficients.get_types();
    std::list<GlossaryType> TestedGlossaryType;
    TestedGlossaryType.assign(vect_types.begin(), vect_types.end());
    expected_types.sort();
    TestedGlossaryType.sort();

    bool expected_types_found = std::ranges::includes(TestedGlossaryType, expected_types);

    MFEM_VERIFY(expected_types_found,
                "Error in " + this->integrator_name_ +
                    " at least one coefficient does not match with the expected list of "
                    "GlossaryType for the current SLOTH-based Integrator. Please check your data.");
  }
}

template <class VARS>
std::optional<Coefficient> SlothNLFormIntegrator<VARS>::get_coefficient(const int blk,
                                                                        GlossaryType type,
                                                                        unsigned int id) {
  Coefficients coefficients = this->coefficients_[blk];

  for (unsigned int i = 0; i < coefficients.size(); i++) {
    auto coef = coefficients[i];
    if (coef.get_type() == type && coef.get_id() == id) return coef;
  }
  return std::nullopt;
}

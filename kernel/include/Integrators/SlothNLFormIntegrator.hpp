/**
 * @file SlothNLFormIntegrator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Class prodiving common methods to all NLFormIntegrators
 * @version 0.1
 * @date 2025-09-05
 *
 * @copyright CEA (C) 2025
 *
 * @anchor SlothNLFormIntegrator
 *
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

#include <algorithm>
#include <list>
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
  std::vector<mfem::ParGridFunction> aux_old_gf_;
  std::vector<std::vector<std::string>> aux_infos_;

  std::vector<VARS*> auxvariables_;
  Parameters params_;
  std::vector<Coefficients> coefficients_;
  unsigned int nb_blk_;

  std::vector<mfem::ParGridFunction> get_aux_gf();
  std::vector<std::vector<std::string>> get_aux_infos();

  void check_coefficient_types(std::list<GlossaryType> expected_type);
  std::optional<Coefficient> get_coefficient(const int blk, GlossaryType type, unsigned int id,
                                             std::optional<int> bdr_id = std::nullopt);

 public:
  virtual void init() = 0;
  SlothNLFormIntegrator(const std::vector<mfem::ParGridFunction> u_old,
                        const std::vector<mfem::ParGridFunction> aux_old, const Parameters& params,
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
                                                   const std::vector<mfem::ParGridFunction> aux_old,
                                                   const Parameters& params,
                                                   std::vector<VARS*> auxvars,
                                                   const std::vector<Coefficients>& coefficients)
    : u_old_(u_old), aux_old_gf_(aux_old), params_(params), coefficients_(coefficients) {
  this->nb_blk_ = this->u_old_.size();
  this->manage_auxiliary_variables(auxvars);
  this->aux_gf_ = this->get_aux_gf();
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
      this->vect_aux_gf_.emplace_back(std::move(auxvar.get_gf()));

      std::vector<std::string> var_info = auxvar.get_additional_variable_info();

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
 * @brief Verify that the coefficients associated with the integrator match the expected glossary
 *        types.
 *
 * The check is order-independent and ignores duplicates. If any required glossary type is missing,
 * the function aborts execution via MFEM_VERIFY.
 *
 * @tparam VARS Template parameter defining the variable set handled by the integrator.
 *
 * @param expected_types List of glossary types required by the current SLOTH-based integrator.
 *
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

/**
 * @brief Retrieve a coefficient by type and identifier.
 *
 * This function searches the list of coefficients associated with the
 * specified block (blk) and returns the coefficient that matches
 * the given glossary type, identifier and if given, the boundary id.
 *
 * @tparam VARS Template parameter defining the variables.
 *
 * @param blk     Index of the block.
 * @param type    Type of coefficient.
 * @param id      Identifier of the coefficient.
 * @param bdr_id  Optional boundary id for boundary coefficients.
 *
 * @return An std::optional containing the matching Coefficient if found;
 *         std::nullopt otherwise.
 */
template <class VARS>
std::optional<Coefficient> SlothNLFormIntegrator<VARS>::get_coefficient(const int blk,
                                                                        GlossaryType type,
                                                                        unsigned int id,
                                                                        std::optional<int> bdr_id) {
  Coefficients coefficients = this->coefficients_[blk];

  for (unsigned int i = 0; i < coefficients.size(); i++) {
    auto coef = coefficients[i];
    if (coef.get_type() == type && coef.get_id() == id) {
      if (bdr_id.has_value()) {
        auto bdr_index = coef.get_bdr_index_coef();
        bool bdr_id_found =
            std::find(bdr_index.begin(), bdr_index.end(), bdr_id) != bdr_index.end();
        if (bdr_id_found) return coef;
      } else {
        return coef;
      }
    }
  }
  return std::nullopt;
}

/**
 * @file MeltingCalphadNLFormIntegrator.tpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief NonlinearFormIntegrator for ad-hoc melting term based on temperature
 * @version 0.1
 * @date 2025-09-05
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
#include <algorithm>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Integrators/MeltingBaseNLFormIntegrator.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new MeltingCalphadNLFormIntegrator object.
 *
 * This constructor initializes the nonlinear form integrator. It forwards the provided previous
 * solution fields, simulation parameters, auxiliary variables, and coefficients to the base melting
 * nonlinear form integrator.
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 *
 * @param u_old        Vector of previous-time-step solution fields.
 * @param params       Paramters that can be used with the integrator.
 * @param auxvars      Auxiliary variables required by the inetgrator.
 * @param coefficients List of coefficients defining material properties.
 *
 */
template <class VARS>
MeltingCalphadNLFormIntegrator<VARS>::MeltingCalphadNLFormIntegrator(
    Geometry geometry, const std::vector<mfem::ParGridFunction>& u_old,
    const std::vector<mfem::ParGridFunction>& aux_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients)
    : MeltingBaseNLFormIntegrator<VARS>(geometry, u_old, aux_old, params, auxvars, coefficients) {
  this->integrator_name_ = "MeltingCalphad";

  this->get_parameters();
  this->check_driving_forces();
  this->check_nucleus();
}

/**
 * @brief Retrieve and initialize model parameters from the parameter set.
 *
 * This function reads the parameters required by the integrator from `params_` and stores them
 * in the corresponding member variables.
 *
 * The following parameters are queried:
 * - `melting_factor` (double, optional):  Scaling factor applied to the melting contribution.
 * Defaults to `1.0` if not provided.
 * - `primary_phase` (std::string, required): Name of the primary phase (eg. SOLID).
 * - `secondary_phase` (std::string, required): Name of the secondary phase (eg. LIQUID).
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 *
 * @pre `params_` must be initialized and contain the required entries.
 */
template <class VARS>
void MeltingCalphadNLFormIntegrator<VARS>::get_parameters() {
  this->alpha_ = this->params_.template get_param_value_or_default<double>("melting_factor", 1.);
  //
  this->primary_phase_ = this->params_.template get_param_value<std::string>("primary_phase");
  this->secondary_phase_ = this->params_.template get_param_value<std::string>("secondary_phase");
}

/**
 * @brief Check the presence of driving forces variables for two phases and stores them in a
 * container 'dgm_'
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 */
template <class VARS>
void MeltingCalphadNLFormIntegrator<VARS>::check_driving_forces() {
  bool primary_phase_found = false;
  bool secondary_phase_found = false;
  this->dgm_primary_phase_index_ = -1;
  this->dgm_secondary_phase_index_ = -1;
  for (std::size_t i = 0; i < this->aux_infos_.size(); ++i) {
    const auto& variable_info = this->aux_infos_[i];
    MFEM_VERIFY(!variable_info.empty(), "Empty variable_info encountered.");

    size_t vsize = variable_info.size();
    MFEM_VERIFY(vsize > 0,
                "MeltingCalphadNLFormIntegrator<VARS, SCHEME, ENERGY, "
                "MOBI,INTERPOLATION>: at least "
                "one additionnal information is expected for auxiliary variables associated with "
                "this integrator");
    const std::string& symbol = toLowerCase(variable_info.back());

    if (calphad_outputs::from(symbol) != calphad_outputs::dgm) continue;
    MFEM_VERIFY(vsize == 2,
                "Error while getting driving forces. Two additional informations are excepted "
                ": the name of the phase and the symbol 'dgm'");

    if (variable_info[0] == this->primary_phase_) {
      this->dgm_primary_phase_index_ = i;
      primary_phase_found = true;
    } else if (variable_info[0] == this->secondary_phase_) {
      this->dgm_secondary_phase_index_ = i;
      secondary_phase_found = true;
    }
  }
  MFEM_VERIFY(primary_phase_found && secondary_phase_found,
              "Both primary and secondary driving forces must be set.");
}

/**
 * @brief Check the presence of nucleus variable  and stores it in a
 * container 'nucleus_'
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 */
template <class VARS>
void MeltingCalphadNLFormIntegrator<VARS>::check_nucleus() {
  bool nucleus_found = false;
  this->nucleus_index_ = -1;
  for (std::size_t i = 0; i < this->aux_infos_.size(); ++i) {
    const auto& variable_info = this->aux_infos_[i];
    MFEM_VERIFY(!variable_info.empty(), "Empty variable_info encountered.");

    size_t vsize = variable_info.size();
    MFEM_VERIFY(vsize > 0,
                "MeltingCalphadNLFormIntegrator<VARS, SCHEME, ENERGY, "
                "MOBI,INTERPOLATION>: at least "
                "one additionnal information is expected for auxiliary variables associated with "
                "this integrator");
    const std::string& symbol = toLowerCase(variable_info.back());

    if (calphad_outputs::from(symbol) != calphad_outputs::nucleus) continue;
    MFEM_VERIFY(vsize == 2,
                "Error while getting nucleus. Two additional informations are excepted "
                ": the name of the phase and the symbol 'nucleus'");

    if (variable_info[0] == this->secondary_phase_) {
      this->nucleus_index_ = i;
      nucleus_found = true;
    }
  }
  MFEM_VERIFY(nucleus_found, "Nucleus for secondary phase must be set.");
}

/**
 * @brief Compute the value of the phase change at integration point
 *
 * @remark Written for two phase
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param Tr Element transformation.
 * @param ir Integration point
 * @param blk   Index of the block.
 * @param u value of the current solution at ip
 * @param un value of the previous solution at ip
 *
 * @return The computed phase change at integration point
 */
template <class VARS>
double MeltingCalphadNLFormIntegrator<VARS>::get_phase_change_at_ip(
    [[maybe_unused]] unsigned int blk, [[maybe_unused]] const std::span<const double>& values,
    [[maybe_unused]] const std::span<const double>& aux_values) {
  double primary_dgm = -1.;
  double secondary_dgm = -1.;
  if (this->dgm_primary_phase_index_ > 0 && this->dgm_secondary_phase_index_ > 0) {
    primary_dgm = aux_values[this->dgm_primary_phase_index_];
    secondary_dgm = aux_values[this->dgm_secondary_phase_index_];
  }
  constexpr double epsilon = 1.e-12;
  return (std::abs(primary_dgm) > epsilon && std::abs(secondary_dgm) > epsilon)
             ? (this->alpha_ * (secondary_dgm - primary_dgm))
             : 0.;
}

/**
 * @brief Compute the value of a seed of the secondary phase at integration point
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param Tr Element transformation.
 * @param ir Integration point
 * @param blk   Index of the block.
 * @param u value of the current solution at ip
 * @param un value of the previous solution at ip
 *
 * @return The computed seed at integration point
 */
template <class VARS>
double MeltingCalphadNLFormIntegrator<VARS>::get_seed_at_ip(
    [[maybe_unused]] unsigned int blk, [[maybe_unused]] const std::span<const double>& values,
    [[maybe_unused]] const std::span<const double>& aux_values) {
  // Nucleus must be equal to zero except when phase transition starts

  // std::span<const double> u(values.begin(), values.begin() + this->nb_blk_);
  // std::span<const double> un(values.begin() + this->nb_blk_, values.end());

  const double seed = -aux_values[this->nucleus_index_];

  return seed;
}

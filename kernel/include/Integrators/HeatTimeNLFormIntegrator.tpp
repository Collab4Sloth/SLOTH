/**
 * @file HeatTimeNLFormIntegrator.tpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief VF of the time derivative
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
#include <string>
#include <utility>
#include <vector>

#include "Coefficients/AxiCylindricalCoefficient.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "Integrators/TimeNLFormIntegrator.hpp"
#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new HeatTimeNLFormIntegrator object.
 *
 * This constructor initializes the nonlinear form integrator.
 * It forwards the provided previous solution fields, simulation
 * parameters, auxiliary variables, and coefficients to the base
 * SLOTH nonlinear form integrator.
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
HeatTimeNLFormIntegrator<VARS>::HeatTimeNLFormIntegrator(
    Geometry geometry, const double time_step, const std::vector<mfem::ParGridFunction> u_old,
    const std::vector<mfem::ParGridFunction>& aux_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients)
    : TimeNLFormIntegrator<VARS>(geometry, time_step, u_old, aux_old, params, auxvars,
                                 coefficients) {
  this->integrator_name_ = "HeatTime";

  this->expected_list_.push_back(GlossaryType::Concentration);
  this->expected_list_.push_back(GlossaryType::HeatCapacity);
}

/**
 * @brief Retrieve and store the coefficients required by the HeatTimeNLFormIntegrator integrator.
 *
 * This method collects two default coefficients:
 * - 'coefficient_A' stores the first coefficient as a Concentration,
 * - 'coefficient_B' stores the second coefficient as a HeatCapacity,
 *
 *
 * @remark By default, A and B are equal to one. This method specializes the coefficients for the
 * heat transfer equation
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 */
template <class VARS>
void HeatTimeNLFormIntegrator<VARS>::get_coefficients() {
  for (unsigned int i = 0; i < this->nb_blk_; i++) {
    if (this->get_coefficient(i, GlossaryType::Concentration, 0).has_value()) {
      this->coefficient_A.add(*(this->get_coefficient(i, GlossaryType::Concentration, 0)));
    }
    if (this->get_coefficient(i, GlossaryType::HeatCapacity, 0).has_value()) {
      this->coefficient_B.add(*(this->get_coefficient(i, GlossaryType::HeatCapacity, 0)));
    }
  }
}

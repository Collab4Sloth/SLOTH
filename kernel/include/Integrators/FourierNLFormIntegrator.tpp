/**
 * @file FourierNLFormIntegrator.tpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief FV for the heat transfer equation
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
#include <vector>

#include "Coefficients/AxiCylindricalCoefficient.hpp"
#include "Integrators/DiffusionNLFormIntegrator.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new FourierNLFormIntegrator object.
 *
 * This constructor initializes a nonlinear form integrator for the
 * heat transfer equation (Fourier's law).
 *
 * It also sets the internal integrator name to "Fourier" and defines
 * the expected list of glossary types required for this integrator.
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 *
 * @param u_old        Vector of previous-time-step solution fields.
 * @param params       Parameters that can be used by the integrator.
 * @param auxvars      Auxiliary variables required for the integrator.
 * @param coefficients List of coefficients defining material
 *                     properties.
 *
 */
template <class VARS>
FourierNLFormIntegrator<VARS>::FourierNLFormIntegrator(
    Geometry geometry, const std::vector<mfem::ParGridFunction>& u_old,
    const std::vector<mfem::ParGridFunction>& aux_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients)
    : DiffusionNLFormIntegrator<VARS>(geometry, u_old, aux_old, params, auxvars, coefficients) {
  this->integrator_name_ = "Fourier";
  this->expected_list_.push_back(GlossaryType::Conductivity);
}

/**
 * @brief Retrieve and store the diffusion coefficients for the integrator.
 *
 * This method collects the coefficients of type `GlossaryType::Conductivity`
 * from each coefficients and adds them to the internal diffusion
 * storage. Only the coefficient with ID 0 is considered for each block.
 *
 * @remark This method can be overridden by derived classes to provide
 *         custom behavior for retrieving coefficients.
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 */
template <class VARS>
void FourierNLFormIntegrator<VARS>::get_coefficients() {
  for (unsigned int i = 0; i < this->nb_blk_; i++) {
    if (this->get_coefficient(i, GlossaryType::Conductivity, 0).has_value()) {
      this->diffusion.add(*(this->get_coefficient(i, GlossaryType::Conductivity, 0)));
    }
  }
}

/**
 * @file NeumannNLFormIntegrator.tpp
 * @author Marine Harel (marine.harel@cea.fr)
 * @brief Neumann integrator
 * @version 0.1
 * @date 2026-02-23
 *
 * @copyright CEA (C) 2026
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

#include <list>
#include <memory>
#include <utility>
#include <vector>

#include "Integrators/NeumannNLFormIntegrator.hpp"
#include "Integrators/RobinNLFormIntegrator.hpp"
#include "Integrators/SlothGridFunction.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameters.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Retrieve and store the coefficients required by the integrator.
 *
 * This method collects a coefficient of type Neumann and stores it.
 *
 * @remark This method can be overridden in derived classes to provide
 *         custom behavior for retrieving coefficients.
 *
 * @tparam VARS Template parameter defining the variabls used in the integrator.
 *
 */
template <class VARS>
void NeumannNLFormIntegrator<VARS>::get_coefficients() {
  auto neumann_opt = this->get_coefficient(this->blk_, GlossaryType::Neumann, 0, this->bdr_id_);
  if (neumann_opt.has_value()) {
    this->robin_b = *neumann_opt;
  }
  // If no neumann coefficient is given, robin_a and robin_b stay zero by default
  // (Homogeneous Neumann)
}

/**
 * @brief Construct a new NeumannNLFormIntegrator object.
 *
 * This constructor initializes the nonlinear form integrator. It forwards the provided previous
 * solution fields, simulation parameters, auxiliary variables, and coefficients to the base SLOTH
 * nonlinear form integrator.
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 *
 * @param u_old        Vector of previous-time-step solution fields.
 * @param params       Paramters that can be used with the integrator.
 * @param auxvars      Auxiliary variables required by the inetgrator.
 * @param coefficients List of coefficients defining material properties.
 * @param block        Block to which the Neumann condition is applied.
 * @param bdr_id       Id of the boundary to which the Neumann condition is applied.
 *
 */
template <class VARS>
NeumannNLFormIntegrator<VARS>::NeumannNLFormIntegrator(
    const std::vector<mfem::ParGridFunction>& u_old,
    const std::vector<mfem::ParGridFunction>& aux_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients,
    const unsigned int block, const unsigned int bdr_id)
    : RobinNLFormIntegrator<VARS>(u_old, aux_old, params, auxvars, coefficients, block, bdr_id) {
  this->integrator_name_ = "Neumann";
}

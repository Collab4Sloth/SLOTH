/**
 * @file AdvectionBoundaryCHNLFormIntegrator.tpp
 * @author Tom Schnetz (tom.schnetz@cea.fr)
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Inflow/outflow advection boundary term for the split Cahn-Hilliard equations
 * @version 0.5
 * @date 2026-10-17
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

#include <string>
#include <utility>
#include <vector>

#include "Integrators/AdvectionBoundaryCHNLFormIntegrator.hpp"
#include "Integrators/AdvectionBoundaryNLFormIntegrator.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new AdvectionBoundaryCHNLFormIntegrator object.
 *
 * Velocity and Inflow of phi (block 0), residual row 1 (equation holding dphi/dt), advected
 * variable 0 (phi).
 *
 * @param block   Block holding the Inflow coefficient (must be phi, block 0).
 * @param bdr_id  Boundary id.
 */
template <class VARS>
AdvectionBoundaryCHNLFormIntegrator<VARS>::AdvectionBoundaryCHNLFormIntegrator(
    Geometry geometry, const double time_step, const std::vector<mfem::ParGridFunction>& u_old,
    const std::vector<mfem::ParGridFunction>& aux_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients,
    const unsigned int block, const unsigned int bdr_id)
    : AdvectionBoundaryNLFormIntegrator<VARS>(geometry, time_step, u_old, aux_old, params,
                                              std::move(auxvars), coefficients, 0, 1, 0, 0,
                                              bdr_id) {
  this->integrator_name_ = "AdvectionBoundaryCH";
  MFEM_VERIFY(block == 0, "AdvectionBoundaryCHNLFormIntegrator: the Inflow coefficient must be "
                          "attached to phi (block 0), not to block " +
                              std::to_string(block) + ".");
}

/**
 * @file AdvectionBoundaryCHNLFormIntegrator.hpp
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

#include <vector>

#include "Integrators/AdvectionBoundaryNLFormIntegrator.hpp"
#include "Parameters/Parameters.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Inflow/outflow boundary term of the advection operator for the split (phi, mu)
 * Cahn-Hilliard formulation: the advection of phi (block 0) is added to the residual of block 1
 * (which holds dphi/dt):
 *   R_1 = ( (v.n)^+ phi + (v.n)^- phi_in, psi_1 )_Gamma
 * The Velocity and Inflow coefficients are those of phi (block 0), consistently with
 * SplitAdvectionNLFormIntegrator.
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 */
template <class VARS>
class AdvectionBoundaryCHNLFormIntegrator : public AdvectionBoundaryNLFormIntegrator<VARS> {
 public:
  AdvectionBoundaryCHNLFormIntegrator(Geometry geometry, const double time_step,
                                      const std::vector<mfem::ParGridFunction>& u_old,
                                      const std::vector<mfem::ParGridFunction>& aux_old,
                                      const Parameters& params, std::vector<VARS*> auxvars,
                                      const std::vector<Coefficients>& coefficients,
                                      const unsigned int block, const unsigned int bdr_id);
};

#include "Integrators/AdvectionBoundaryCHNLFormIntegrator.tpp"

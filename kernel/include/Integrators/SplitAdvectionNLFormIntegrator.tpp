/**
 * @file SplitAdvectionNLFormIntegrator.tpp
 * @author Tom Schnetz (tom.schnetz@cea.fr)
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief VF of the advection term in the split Cahn-Hilliard equations
 * @version 0.5
 * @date 2026-10-17
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
#include <utility>
#include <vector>

#include "Integrators/AdvectionNLFormIntegrator.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Construct a new SplitAdvectionNLFormIntegrator object (same arguments as Advection).
 */
template <class VARS>
SplitAdvectionNLFormIntegrator<VARS>::SplitAdvectionNLFormIntegrator(
    Geometry geometry, const double time_step, const std::vector<mfem::ParGridFunction>& u_old,
    const std::vector<mfem::ParGridFunction>& aux_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients)
    : AdvectionNLFormIntegrator<VARS>(geometry, time_step, u_old, aux_old, params,
                                      std::move(auxvars), coefficients) {
  this->integrator_name_ = "SplitAdvection";
}

/**
 * @brief Single term: residual of block 1, advected variable phi (block 0), velocity of phi.
 */
template <class VARS>
void SplitAdvectionNLFormIntegrator<VARS>::build_terms() {
  MFEM_VERIFY(this->nb_blk_ == 2,
              "SplitAdvectionNLFormIntegrator: two blocks (phi, mu) are expected.");
  MFEM_VERIFY(this->velocity_[0].active(),
              "SplitAdvectionNLFormIntegrator: Velocity (or VelocityPrefactor) coefficients are "
              "required for phi (block 0).");
  this->terms_.clear();
  this->terms_.push_back({1, 0, 0});
}

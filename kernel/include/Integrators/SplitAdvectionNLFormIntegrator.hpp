/**
 * @file SplitAdvectionNLFormIntegrator.hpp
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
#include <vector>

#include "Integrators/AdvectionNLFormIntegrator.hpp"
#include "Parameters/Parameters.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Advection term for the split (phi, mu) Cahn-Hilliard formulation.
 *
 * Single term (r, c) = (1, 0): the advection of phi (block 0) is added to the residual of
 * block 1, which holds dphi/dt (SplitTimeDerivative):
 *   R_1 = - ( phi v, grad(psi_1) ).
 * The velocity is given by the Velocity coefficients of block 0 (phi).
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 */
template <class VARS>
class SplitAdvectionNLFormIntegrator : public AdvectionNLFormIntegrator<VARS> {
 protected:
  void build_terms() override;

 public:
  SplitAdvectionNLFormIntegrator(Geometry geometry, const double time_step,
                                 const std::vector<mfem::ParGridFunction>& u_old,
                                 const std::vector<mfem::ParGridFunction>& aux_old,
                                 const Parameters& params, std::vector<VARS*> auxvars,
                                 const std::vector<Coefficients>& coefficients);
};

#include "Integrators/SplitAdvectionNLFormIntegrator.tpp"

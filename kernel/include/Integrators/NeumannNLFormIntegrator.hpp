/**
 * @file NeumannNLFormIntegrator.hpp
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

#include "Integrators/RobinNLFormIntegrator.hpp"
#include "Integrators/SlothGridFunction.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameters.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief  Class dedicated to the Neumann boundary condiiton integrator
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 */
template <class VARS>
class NeumannNLFormIntegrator : public RobinNLFormIntegrator<VARS> {
 protected:
  void get_coefficients() override;

 public:
  NeumannNLFormIntegrator(Geometry geometry, const std::vector<mfem::ParGridFunction>& u_old,
                          const std::vector<mfem::ParGridFunction>& aux_old,
                          const Parameters& params, std::vector<VARS*> auxvars,
                          const std::vector<Coefficients>& coefficients, const unsigned int block,
                          const unsigned int bdr_id);

  virtual ~NeumannNLFormIntegrator() = default;
};

#include "Integrators/NeumannNLFormIntegrator.tpp"

/**
 * @file FickNLFormIntegrator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief FV for the mass diffusion equation
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

#include "Integrators/DiffusionNLFormIntegrator.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @class FickNLFormIntegrator
 * @brief Nonlinear form integrator for the right-hand side of the mass diffusion equation
 *        (Fick's law).
 *
 * This integrator assembles the nonlinear form corresponding to the
 * diffusive mass flux governed by Fick's law. It extends the
 * DiffusionNLFormIntegrator by providing the appropriate coefficients
 * required for mass diffusion problems.
 *
 * @tparam VARS Template parameter defining the variables used by the
 *              integrator.
 */
template <class VARS>
class FickNLFormIntegrator : public DiffusionNLFormIntegrator<VARS> {
 protected:
  void get_coefficients() override;

 public:
  FickNLFormIntegrator(const std::vector<mfem::ParGridFunction> u_old,
                       const std::vector<mfem::ParGridFunction>& aux_old, const Parameters& params,
                       std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients);
  virtual ~FickNLFormIntegrator() = default;
};

#include "Integrators/FickNLFormIntegrator.tpp"

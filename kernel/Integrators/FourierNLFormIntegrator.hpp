/**
 * @file FourierNLFormIntegrator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief FV for the heat transfer equation
 * @version 0.1
 * @date 2025-09-05
 *
 * Copyright CEA (C) 2025
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

#pragma once

/**
 * @brief  Class dedicated to the FV of the mass diffusion equation
 *
 * @tparam SCHEME
 * @tparam DIFFU_NAME
 */
template <class VARS>
class FourierNLFormIntegrator : public DiffusionNLFormIntegrator<VARS> {
 protected:
  void get_coefficients() override;

 public:
  FourierNLFormIntegrator(const std::vector<mfem::ParGridFunction>& u_old, const Parameters& params,
                          std::vector<VARS*> auxvars,
                          const std::vector<Coefficients>& coefficients);
  virtual ~FourierNLFormIntegrator() = default;
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Construct a new FourierNLFormIntegrator<SCHEME, COEFFICIENT>::FourierNLFormIntegrator
 * object
 *
 * @tparam SCHEME
 * @tparam COEFFICIENT
 * @param u_old
 * @param alpha
 * @param kappa
 */
template <class VARS>
FourierNLFormIntegrator<VARS>::FourierNLFormIntegrator(
    const std::vector<mfem::ParGridFunction>& u_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients)
    : DiffusionNLFormIntegrator<VARS>(u_old, params, auxvars, coefficients) {
  this->integrator_name_ = "Fourier";
  this->expected_list_.push_back(GlossaryType::Conductivity);
}

/**
 * @brief Get diffusion coefficients
 * @remark Could be overridden by child classes
 *
 * @tparam VARS
 */
template <class VARS>
void FourierNLFormIntegrator<VARS>::get_coefficients() {
  for (unsigned int i = 0; i < this->nb_blk_; i++) {
    if (this->get_coefficient(i, GlossaryType::Conductivity, 0).has_value()) {
      this->diffusion.add(*(this->get_coefficient(i, GlossaryType::Conductivity, 0)));
    }
  }
}
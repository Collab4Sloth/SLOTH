/**
 * @file HeatTimeNLFormIntegrator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief VF of the time derivative
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
#include <string>
#include <utility>
#include <vector>

#include "Integrators/SlothNLFormIntegrator.hpp"
#include "Integrators/TimeNLFormIntegrator.hpp"
#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

/**
 * @brief Class dedicated to the FV of the Allen Cahn equation
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 */
template <class VARS>
class HeatTimeNLFormIntegrator : public TimeNLFormIntegrator<VARS> {
 protected:
  void get_coefficients() override;

 public:
  HeatTimeNLFormIntegrator(const std::vector<mfem::ParGridFunction> u_old, const Parameters& params,
                           std::vector<VARS*> auxvars,
                           const std::vector<Coefficients>& coefficients);
  virtual ~HeatTimeNLFormIntegrator() = default;
};
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Construct a new HeatTimeNLFormIntegrator object
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @param u_old
 * @param omega
 * @param lambda
 * @param mob
 */
template <class VARS>
HeatTimeNLFormIntegrator<VARS>::HeatTimeNLFormIntegrator(
    const std::vector<mfem::ParGridFunction> u_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients)
    : TimeNLFormIntegrator<VARS>(u_old, params, auxvars, coefficients) {
  this->expected_list_.push_back(GlossaryType::Concentration);
  this->expected_list_.push_back(GlossaryType::HeatCapacity);
}

/**
 * @brief Get default coefficients
 * @remark Could be overridden by child classes
 *
 * @tparam VARS
 */
template <class VARS>
void HeatTimeNLFormIntegrator<VARS>::get_coefficients() {
  for (unsigned int i = 0; i < this->nb_blk_; i++) {
    if (this->get_coefficient(i, GlossaryType::Conductivity, 0).has_value()) {
      this->coefficient_A.add(*(this->get_coefficient(i, GlossaryType::Concentration, 0)));
      this->coefficient_B.add(*(this->get_coefficient(i, GlossaryType::HeatCapacity, 0)));
    }
  }
}

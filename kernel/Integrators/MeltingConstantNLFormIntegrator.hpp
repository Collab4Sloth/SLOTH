/**
 * @file MeltingConstantNLFormIntegrator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief NonlinearFormIntegrator for ad-hoc constant melting term
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
#include <utility>
#include <vector>

#include "Integrators/SlothNLFormIntegrator.hpp"
#include "MAToolsProfiling/MATimersAPI.hxx"
#include "MeltingBaseNLFormIntegrator.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

/**
 * @brief NonlinearFormIntegrator for ad-hoc melting term based on temperature
 *
 * @tparam VARS
 */
template <class VARS>
class MeltingConstantNLFormIntegrator : public MeltingBaseNLFormIntegrator<VARS> {
 private:
  double alpha_;

  void get_parameters();

 protected:
  double get_phase_change_at_ip(mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ir,
                                unsigned int blk, const double u, const double un) override;

 public:
  MeltingConstantNLFormIntegrator(const std::vector<mfem::ParGridFunction>& u_old,
                                  const Parameters& params, std::vector<VARS*> auxvars,
                                  const std::vector<Coefficients>& coefficients);

  virtual ~MeltingConstantNLFormIntegrator() = default;
};
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Construct a new Melting Temperature N L Form Integrator< V A R S>:: Melting Temperature N
 * L Form Integrator object
 *
 * @tparam VARS
 * @param u_old
 * @param params
 * @param auxvars
 * @param coefficients
 */
template <class VARS>
MeltingConstantNLFormIntegrator<VARS>::MeltingConstantNLFormIntegrator(
    const std::vector<mfem::ParGridFunction>& u_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients)
    : MeltingBaseNLFormIntegrator<VARS>(u_old, params, auxvars, coefficients) {
  this->integrator_name_ = "MeltingConstant";

  this->get_parameters();
}

/**
 * @brief Get parameters
 *
 * @tparam VARS
 */
template <class VARS>
void MeltingConstantNLFormIntegrator<VARS>::get_parameters() {
  this->alpha_ = this->params_.template get_param_value_or_default<double>("melting_factor", 1.);
}

/**
 * @brief Get the value of the phae change at integration point
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam INTERPOLATION
 * @param Tr
 * @param ir
 * @return double
 */
template <class VARS>
double MeltingConstantNLFormIntegrator<VARS>::get_phase_change_at_ip(
    [[maybe_unused]] mfem::ElementTransformation& Tr,
    [[maybe_unused]] const mfem::IntegrationPoint& ir, [[maybe_unused]] unsigned int blk,
    [[maybe_unused]] const double u, [[maybe_unused]] const double un) {
  return this->alpha_;
}

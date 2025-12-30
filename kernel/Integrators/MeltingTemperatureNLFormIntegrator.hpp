/**
 * @file MeltingTemperatureNLFormIntegrator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief NonlinearFormIntegrator for ad-hoc melting term based on temperature
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
class MeltingTemperatureNLFormIntegrator : public MeltingBaseNLFormIntegrator<VARS> {
 private:
  std::vector<mfem::ParGridFunction> temp_gf_;

  double melting_temperature_;
  double melting_enthalpy_;
  void get_parameters();
  void check_variables_consistency();

 protected:
  double get_phase_change_at_ip(mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ir,
                                unsigned int blk, const double u, const double un) override;

 public:
  MeltingTemperatureNLFormIntegrator(const std::vector<mfem::ParGridFunction>& u_old,
                                     const Parameters& params, std::vector<VARS*> auxvars,
                                     const std::vector<Coefficients>& coefficients);

  virtual ~MeltingTemperatureNLFormIntegrator() = default;
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
MeltingTemperatureNLFormIntegrator<VARS>::MeltingTemperatureNLFormIntegrator(
    const std::vector<mfem::ParGridFunction>& u_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients)
    : MeltingBaseNLFormIntegrator<VARS>(u_old, params, auxvars, coefficients) {
  this->integrator_name_ = "MeltingTemperature";
  this->get_parameters();
  this->check_variables_consistency();
}

/**
 * @brief Check variables consistency
 *
 * @tparam VARS
 */
template <class VARS>
void MeltingTemperatureNLFormIntegrator<VARS>::check_variables_consistency() {
  // Temperature scaling for mobility
  for (std::size_t i = 0; i < this->aux_infos_.size(); ++i) {
    const auto& variable_info = this->aux_infos_[i];
    MFEM_VERIFY(!variable_info.empty(), "Empty variable_info encountered.");
    size_t vsize = variable_info.size();

    MFEM_VERIFY(vsize >= 1,
                "MeltingTemperatureNLFormIntegrator<VARS>: at least "
                "one additionnal information is expected for auxiliary variables associated with "
                "this integrator");
    const std::string& symbol = toUpperCase(variable_info.back());
    if (symbol == "T") {
      this->temp_gf_.emplace_back(std::move(this->aux_gf_[i]));
      break;
    }
  }
}

/**
 * @brief Get parameters
 *
 * @tparam VARS
 */
template <class VARS>
void MeltingTemperatureNLFormIntegrator<VARS>::get_parameters() {
  this->melting_temperature_ =
      this->params_.template get_param_value<double>("melting_temperature");
  this->melting_enthalpy_ = this->params_.template get_param_value<double>("melting_enthalpy");
}

/**
 * @brief Compute phase change at integration point
 *
 * @tparam VARS
 * @param Tr
 * @param ir
 * @param blk
 * @param u
 * @param un
 * @return double
 */
template <class VARS>
double MeltingTemperatureNLFormIntegrator<VARS>::get_phase_change_at_ip(
    mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ir,
    [[maybe_unused]] unsigned int blk, [[maybe_unused]] const double u,
    [[maybe_unused]] const double un) {
  const double temperature_at_ip = this->temp_gf_[0].GetValue(Tr, ir);
  double phase_change_at_ip = 0.;
  if (temperature_at_ip > this->melting_temperature_) {
    phase_change_at_ip = this->melting_enthalpy_;
  }
  return phase_change_at_ip;
}

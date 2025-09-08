/**
 * @file AllenCahnTemperatureMeltingNLFormIntegrator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief FV for the Allen Cahn equation with ad-hoc melting controled by temperature and a given
 * melting temperature
 * @todo Extend coefficient to melting temperature and enthalpy
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
#include <memory>
#include <string>
#include <vector>

#include "Integrators/AllenCahnMeltingBaseNLFormIntegrator.hpp"
#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
class AllenCahnTemperatureMeltingNLFormIntegrator final
    : public AllenCahnMeltingBaseNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI, INTERPOLATION> {
 private:
  double melting_temperature_;
  double melting_enthalpy_;
  void get_parameters();
  void check_variables_consistency();

 protected:
  double get_phase_change_at_ip(mfem::ElementTransformation& Tr,
                                const mfem::IntegrationPoint& ir) override;

 public:
  AllenCahnTemperatureMeltingNLFormIntegrator(const std::vector<mfem::ParGridFunction>& u_old,
                                              const Parameters& params, std::vector<VARS*> auxvars);
  ~AllenCahnTemperatureMeltingNLFormIntegrator();
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Construct AllenCahnTemperatureMeltingNLFormIntegrator object
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam INTERPOLATION
 * @param u_old
 * @param params
 * @param aux_gf
 * @return AllenCahnTemperatureMeltingNLFormIntegrator<VARS,SCHEME, ENERGY, MOBI, INTERPOLATION>::
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
AllenCahnTemperatureMeltingNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI, INTERPOLATION>::
    AllenCahnTemperatureMeltingNLFormIntegrator(const std::vector<mfem::ParGridFunction>& u_old,
                                                const Parameters& params,
                                                std::vector<VARS*> auxvars)
    : AllenCahnMeltingBaseNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI, INTERPOLATION>(u_old, params,
                                                                                      auxvars) {
  this->get_parameters();
  this->check_variables_consistency();
}

/**
 * @brief Check variables consistency
 *
 * @tparam VARS
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam INTERPOLATION
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
void AllenCahnTemperatureMeltingNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI,
                                                 INTERPOLATION>::check_variables_consistency() {
  this->aux_gf_infos_ = this->get_aux_infos();

  bool temperature_found = false;
  for (std::size_t i = 0; i < this->aux_gf_infos_.size(); ++i) {
    const auto& variable_info = this->aux_gf_infos_[i];
    MFEM_VERIFY(!variable_info.empty(), "Empty variable_info encountered.");
    size_t vsize = variable_info.size();

    MFEM_VERIFY(vsize >= 1,
                "AllenCahnTemperatureMeltingNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI, "
                "INTERPOLATION>: at least "
                "one additionnal information is expected for auxiliary variables associated with "
                "this integrator");
    const std::string& symbol = toUpperCase(variable_info.back());
    if (symbol == "T") {
      temperature_found = true;
      break;
    }
  }
  MFEM_VERIFY(temperature_found,
              "AllenCahnTemperatureMeltingNLFormIntegrator: "
              "Temperature variable is required as auxiliary variable for this integrator");
}

/**
 * @brief Get parameters
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam INTERPOLATION
 * @param params
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
void AllenCahnTemperatureMeltingNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI,
                                                 INTERPOLATION>::get_parameters() {
  this->melting_temperature_ =
      this->params_.template get_param_value<double>("melting_temperature");
  this->melting_enthalpy_ = this->params_.template get_param_value<double>("melting_enthalpy");
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
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
double AllenCahnTemperatureMeltingNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI, INTERPOLATION>::
    get_phase_change_at_ip(mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ir) {
  const double temperature_at_ip = this->temp_gf_[0].GetValue(Tr, ir);
  double phase_change_at_ip = 0.;
  if (temperature_at_ip > this->melting_temperature_) {
    phase_change_at_ip = this->melting_enthalpy_;
  }
  return phase_change_at_ip;
}
/**
 * @brief Destroy the AAllenCahnTemperatureMeltingNLFormIntegrator object
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam INTERPOLATION
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
AllenCahnTemperatureMeltingNLFormIntegrator<
    VARS, SCHEME, ENERGY, MOBI, INTERPOLATION>::~AllenCahnTemperatureMeltingNLFormIntegrator() {}

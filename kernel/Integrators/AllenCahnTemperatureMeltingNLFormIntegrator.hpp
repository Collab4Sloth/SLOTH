/**
 * @file AllenCahnTemperatureMeltingNLFormIntegrator.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief FV for the Allen Cahn equation with ad-hoc melting controled by temperature and a given
 * melting temperature
 * @todo Extend coefficient to melting temperature and enthalpy
 * @version 0.1
 * @date 2024-9-3
 *
 * Copyright CEA (c) 2024
 *
 */
#include <algorithm>
#include <tuple>
#include <vector>

#include "Integrators/AllenCahnMeltingBaseNLFormIntegrator.hpp"
#include "Profiling/Profiling.hpp"
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

 protected:
  void get_parameters(const Parameters& vectr_param) override;
  double get_phase_change_at_ip(mfem::ElementTransformation& Tr,
                                const mfem::IntegrationPoint& ir) override;

 public:
  AllenCahnTemperatureMeltingNLFormIntegrator(const mfem::ParGridFunction& _u_old,
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
    AllenCahnTemperatureMeltingNLFormIntegrator(const mfem::ParGridFunction& u_old,
                                                const Parameters& params,
                                                std::vector<VARS*> auxvars)
    : AllenCahnMeltingBaseNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI, INTERPOLATION>(u_old, params,
                                                                                      auxvars) {
  MFEM_VERIFY(this->aux_gf_.size() == 1,
              "AllenCahnTemperatureMeltingNLFormIntegrator requires temperature as the only "
              "auxiliary variable.");

  this->get_parameters(params);
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
void AllenCahnTemperatureMeltingNLFormIntegrator<
    VARS, SCHEME, ENERGY, MOBI, INTERPOLATION>::get_parameters(const Parameters& params) {
  AllenCahnMeltingBaseNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI, INTERPOLATION>::get_parameters(
      params);
  this->melting_temperature_ = params.get_param_value<double>("melting_temperature");
  this->melting_enthalpy_ = params.get_param_value<double>("melting_enthalpy");
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
  const double temperature_at_ip = this->aux_gf_[0].GetValue(Tr, ir);
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

/**
 * @file AllenCahnCalphadMeltingNLFormIntegrator.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Compute the enthalpy of phase change based on two driving forces defined as auxiliary
 * variables
 * @version 0.1
 * @date 2025-05-07
 *
 * Copyright CEA (c) 2025
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
class AllenCahnCalphadMeltingNLFormIntegrator final
    : public AllenCahnMeltingBaseNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI, INTERPOLATION> {
 private:
  double alpha_;
  std::string primary_phase_;
  std::string secondary_phase_;
  std::vector<mfem::Vector> dgm_;
  void get_parameters();

  void check_driving_forces();

 protected:
  double get_phase_change_at_ip(mfem::ElementTransformation& Tr,
                                const mfem::IntegrationPoint& ir) override;

 public:
  AllenCahnCalphadMeltingNLFormIntegrator(const mfem::ParGridFunction& _u_old,
                                          const Parameters& params, std::vector<VARS*> auxvars);
  ~AllenCahnCalphadMeltingNLFormIntegrator();
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Construct AllenCahnCalphadMeltingNLFormIntegrator object
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam INTERPOLATION
 * @param u_old
 * @param params
 * @param aux_gf
 * @return AllenCahnCalphadMeltingNLFormIntegrator<VARS,SCHEME, ENERGY, MOBI, INTERPOLATION>::
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
AllenCahnCalphadMeltingNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI, INTERPOLATION>::
    AllenCahnCalphadMeltingNLFormIntegrator(const mfem::ParGridFunction& u_old,
                                            const Parameters& params, std::vector<VARS*> auxvars)
    : AllenCahnMeltingBaseNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI, INTERPOLATION>(u_old, params,
                                                                                      auxvars) {
  this->get_parameters();
  this->check_driving_forces();
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
void AllenCahnCalphadMeltingNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI,
                                             INTERPOLATION>::get_parameters() {
  this->alpha_ = this->params_.template get_param_value<double>("melting_factor");
  //
  this->primary_phase_ = this->params_.template get_param_value<std::string>("primary_phase");
  this->secondary_phase_ = this->params_.template get_param_value<std::string>("secondary_phase");
}

/**
 * @brief Check presence of driving force
 *
 * @tparam VARS
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam INTERPOLATION
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
void AllenCahnCalphadMeltingNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI,
                                             INTERPOLATION>::check_driving_forces() {
  bool primary_phase_found = false;
  bool secondary_phase_found = false;

  this->dgm_.resize(2);

  for (std::size_t i = 0; i < this->aux_gf_infos_.size(); ++i) {
    const auto& variable_info = this->aux_gf_infos_[i];
    const std::string& symbol = toLowerCase(variable_info.back());

    if (calphad_outputs::from(symbol) != calphad_outputs::dgm) continue;
    // TODO(cci) maybe it could be better to have a variable already calcukated for Deltadgm
    MFEM_VERIFY(variable_info.size() == 2,
                "Error while getting driving forces. Two additional informations are excepted "
                ": the name of the phase and the symbol 'dgm'");

    if (variable_info[0] == this->primary_phase_) {
      this->dgm_[0] = this->aux_gf[i];
      primary_phase_found = true;
    } else if (variable_info[0] == this->secondary_phase_) {
      this->dgm_[1] = this->aux_gf[i];
      secondary_phase_found = true;
    }
  }
  MFEM_VERIFY(primary_phase_found && secondary_phase_found,
              "Both primary and secondary driving forces must be set.");
}
/**
 * @brief Get the value of the phae change at integration point
 * @remark Written for two phase
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
double AllenCahnCalphadMeltingNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI, INTERPOLATION>::
    get_phase_change_at_ip(mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ir) {
  double phase_change_at_ip = this->dgm_[0].GetValue(Tr, ir) - this->dgm_[1].GetValue(Tr, ir);
  return phase_change_at_ip;
}
/**
 * @brief Destroy the AAllenCahnCalphadMeltingNLFormIntegrator object
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam INTERPOLATION
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
AllenCahnCalphadMeltingNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI,
                                        INTERPOLATION>::~AllenCahnCalphadMeltingNLFormIntegrator() {
}

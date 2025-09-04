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
#include <string>
#include <vector>

#include "Integrators/AllenCahnMeltingBaseNLFormIntegrator.hpp"
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
  std::vector<mfem::ParGridFunction> dgm_;
  std::vector<mfem::ParGridFunction> nucleus_;
  void get_parameters();

  void check_driving_forces();

  void check_nucleus();

 protected:
  double get_phase_change_at_ip(mfem::ElementTransformation& Tr,
                                const mfem::IntegrationPoint& ir) override;
  double get_seed_at_ip(mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ir) override;

 public:
  AllenCahnCalphadMeltingNLFormIntegrator(const std::vector<mfem::ParGridFunction>& u_old,
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
    AllenCahnCalphadMeltingNLFormIntegrator(const std::vector<mfem::ParGridFunction>& u_old,
                                            const Parameters& params, std::vector<VARS*> auxvars)
    : AllenCahnMeltingBaseNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI, INTERPOLATION>(u_old, params,
                                                                                      auxvars) {
  this->get_parameters();

  this->check_driving_forces();

  this->check_nucleus();
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
  this->alpha_ = this->params_.template get_param_value_or_default<double>("melting_factor", 1.);
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

  for (std::size_t i = 0; i < this->aux_gf_infos_.size(); ++i) {
    const auto& variable_info = this->aux_gf_infos_[i];
    MFEM_VERIFY(!variable_info.empty(), "Empty variable_info encountered.");

    size_t vsize = variable_info.size();
    MFEM_VERIFY(vsize > 0,
                "AllenCahnCalphadMeltingNLFormIntegrator<VARS, SCHEME, ENERGY, "
                "MOBI,INTERPOLATION>: at least "
                "one additionnal information is expected for auxiliary variables associated with "
                "this integrator");
    const std::string& symbol = toLowerCase(variable_info.back());

    if (calphad_outputs::from(symbol) != calphad_outputs::dgm) continue;
    MFEM_VERIFY(vsize == 2,
                "Error while getting driving forces. Two additional informations are excepted "
                ": the name of the phase and the symbol 'dgm'");

    if (variable_info[0] == this->primary_phase_) {
      this->dgm_.insert(this->dgm_.begin(), this->aux_gf_[i]);
      primary_phase_found = true;
    } else if (variable_info[0] == this->secondary_phase_) {
      this->dgm_.emplace_back(this->aux_gf_[i]);
      secondary_phase_found = true;
    }
  }
  MFEM_VERIFY(primary_phase_found && secondary_phase_found,
              "Both primary and secondary driving forces must be set.");
}

/**
 * @brief Check presence of nucleus
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
                                             INTERPOLATION>::check_nucleus() {
  bool nucleus_found = false;

  for (std::size_t i = 0; i < this->aux_gf_infos_.size(); ++i) {
    const auto& variable_info = this->aux_gf_infos_[i];
    MFEM_VERIFY(!variable_info.empty(), "Empty variable_info encountered.");

    size_t vsize = variable_info.size();
    MFEM_VERIFY(vsize > 0,
                "AllenCahnCalphadMeltingNLFormIntegrator<VARS, SCHEME, ENERGY, "
                "MOBI,INTERPOLATION>: at least "
                "one additionnal information is expected for auxiliary variables associated with "
                "this integrator");
    const std::string& symbol = toLowerCase(variable_info.back());

    if (calphad_outputs::from(symbol) != calphad_outputs::nucleus) continue;
    MFEM_VERIFY(vsize == 2,
                "Error while getting nucleus. Two additional informations are excepted "
                ": the name of the phase and the symbol 'nucleus'");

    if (variable_info[0] == this->secondary_phase_) {
      this->nucleus_.emplace_back(this->aux_gf_[i]);
      nucleus_found = true;
    }
  }
  MFEM_VERIFY(nucleus_found, "Nucleus for secondary phase must be set.");
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
  double primary_dgm = -1.;
  double secondary_dgm = -1.;
  if (this->dgm_.size() == 2) {
    primary_dgm = this->dgm_[0].GetValue(Tr, ir);
    secondary_dgm = this->dgm_[1].GetValue(Tr, ir);
  }
  constexpr double epsilon = 1.e-12;
  return (std::abs(primary_dgm) > epsilon && std::abs(secondary_dgm) > epsilon)
             ? (this->alpha_ * (secondary_dgm - primary_dgm))
             : 0.;
}

template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI, ThermodynamicsPotentials INTERPOLATION>
double
AllenCahnCalphadMeltingNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI, INTERPOLATION>::get_seed_at_ip(
    mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ir) {
  // Nucleus must be equal to zero except when phase transition starts
  const double seed = -this->nucleus_[0].GetValue(Tr, ir);

  return seed;
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

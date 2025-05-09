/**
 * @file MassDiffusionFluxNLFormIntegrator.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Class dedicated to the VF of  diffusion flux (gradient of chemical potential) used in mass
 * balance equation
 * @version 0.1
 * @date 2025-04-03
 *
 * @copyright Copyright (c) 2025
 *
 */

#include <algorithm>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Coefficients/DiffusionCoefficient.hpp"
#include "Coefficients/PhaseFieldPotentials.hpp"
#include "Integrators/SlothGridFunction.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

/**
 * @brief  Class dedicated to the VF of  diffusion flux (gradient of chemical potential) used in
 * mass balance equation
 *
 * @tparam VARS
 */
template <class VARS>
class MassDiffusionFluxNLFormIntegrator : public DiffusionFluxNLFormIntegrator<VARS> {
 private:
  void check_variables_consistency();

 protected:
  std::map<std::string, mfem::ParGridFunction> mu_gf_;
  std::map<std::string, mfem::ParGridFunction> dmu_gf_;
  std::map<std::string, mfem::ParGridFunction> mob_gf_;
  mfem::ParGridFunction temp_gf_;
  bool scaling_by_temperature_{false};

  void get_parameters() override;
  std::vector<mfem::Vector> get_flux_gradient() override;
  std::vector<mfem::Coefficient> get_flux_coefficient() override;

 public:
  MassDiffusionFluxNLFormIntegrator(const mfem::ParGridFunction& u_old, const Parameters& params,
                                    std::vector<VARS*> auxvars);
  ~MassDiffusionFluxNLFormIntegrator();
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Return parameters required by these integrators
 *
 * @tparam VARS
 */
template <class VARS>
void MassDiffusionFluxNLFormIntegrator<VARS>::get_parameters() {
  DiffusionFluxNLFormIntegrator<VARS>::get_parameters();

  // Scaling interdiffusion coefficient by RT
  this->scaling_by_temperature_ =
      this->params_.template get_param_value_or_default<bool>("ScalingByTemperature", false);
}

/**
 * @brief Construct a new MassDiffusionFluxNLFormIntegrator<VARS>::MassDiffusionFluxNLFormIntegrator
 * object
 *
 * @tparam VARS
 * @param u_old
 * @param params
 * @param auxvars
 */
template <class VARS>
MassDiffusionFluxNLFormIntegrator<VARS>::MassDiffusionFluxNLFormIntegrator(
    const mfem::ParGridFunction& u_old, const Parameters& params, std::vector<VARS*> auxvars)
    : DiffusionFluxNLFormIntegrator<VARS>(u_old, params, auxvars) {
  this->check_variables_consistency();
}

/**
 * @brief Check variables consistency
 *
 * @tparam VARS
 */
template <class VARS>
void MassDiffusionFluxNLFormIntegrator<VARS>::check_variables_consistency() {
  std::vector<mfem::ParGridFunction> aux_gf = this->get_aux_gf();
  std::vector<std::vector<std::string>> aux_infos = this->get_aux_infos();

  //==========================================================
  // Get Chemical potentials and mobilities (aux. variables)
  //==========================================================
  bool temperature_found = false;

  for (std::size_t i = 0; i < aux_infos.size(); ++i) {
    const auto& variable_info = aux_infos[i];

    MFEM_VERIFY(!variable_info.empty(), "Empty variable_info encountered.");
    const std::string& symbol = toLowerCase(variable_info.back());

    if (symbol == "mu") {
      MFEM_VERIFY(
          variable_info.size() == 2,
          "MassDiffusionFluxNLFormIntegrator<VARS>::check_variables_consistency: error while "
          "getting chemical potentials. Expected [element, 'mu']");

      const std::string& elem_name = variable_info[0];
      this->mu_gf_.emplace(toUpperCase(elem_name), std::move(aux_gf[i]));
    } else if (symbol == "dmu") {
      MFEM_VERIFY(
          variable_info.size() == 2,
          "MassDiffusionFluxNLFormIntegrator<VARS>::check_variables_consistency: error while "
          "getting diffusion chemical potentials. Expected [element, 'dmu']");

      const std::string& elem_name = variable_info[0];
      this->dmu_gf_.emplace(toUpperCase(elem_name), std::move(aux_gf[i]));
    } else if (symbol == "mob") {
      MFEM_VERIFY(
          variable_info.size() > 2,
          "MassDiffusionFluxNLFormIntegrator<VARS>::check_variables_consistency: error while "
          "getting mobilities.Expected at least [element, 'mob']");

      const std::string& elem_name = variable_info[variable_info.size() - 2];
      this->mob_gf_.emplace(toUpperCase(elem_name), std::move(aux_gf[i]));
    } else if (symbol == "Temperature") {
      this->temp_gf_ = aux_gf[i];
      temperature_found = true;
    }
  }

  MFEM_VERIFY(
      !this->mu_gf_.empty() || !this->dmu_gf_.empty(),
      "MassDiffusionFluxNLFormIntegrator<VARS>::check_variables_consistency: "
      "Neither chemical potentials, nor diffusion chemical potentials found. At least one of them "
      "is required.");

  MFEM_VERIFY(
      this->mu_gf_.empty() || this->dmu_gf_.empty(),
      "MassDiffusionFluxNLFormIntegrator<VARS>::check_variables_consistency: "
      "Either chemical potentials or diffusion chemical potentials can be used, but not both");

  MFEM_VERIFY(
      !this->mob_gf_.empty(),
      "MassDiffusionFluxNLFormIntegrator<VARS>::check_variables_consistency: "
      "No mobilities found. At least as many mobilities as chemical potentials are expected.");

  MFEM_VERIFY(!temperature_found && this->scaling_by_temperature_,
              "MassDiffusionFluxNLFormIntegrator<VARS>::check_variables_consistency: "
              "No temperature found. Expected temperature when ScalingByTemperature is required. ");
}

/**
 * @brief Destroy the MassDiffusionFluxNLFormIntegrator< VARS>::MassDiffusionFluxNLFormIntegrator
 * object
 *
 * @tparam VARS
 */
template <class VARS>
MassDiffusionFluxNLFormIntegrator<VARS>::~MassDiffusionFluxNLFormIntegrator() {}

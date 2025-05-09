/**
 * @file ThermalDiffusionFluxNLFormIntegrator.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief  Class dedicated to the VF of thermal diffusion flux used in energy balance equation
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
 * @brief  Class dedicated to the VF of thermal diffusion flux used in energy balance equation
 *
 * @tparam VARS
 */
template <class VARS>
class ThermalDiffusionFluxNLFormIntegrator : public DiffusionFluxNLFormIntegrator<VARS> {
 private:
  void check_variables_consistency();

 protected:
  mfem::ParGridFunction temp_gf_;
  mfem::ParGridFunction diffu_gf_;

  void get_parameters() override;
  std::vector<mfem::Vector> get_flux_gradient(mfem::ElementTransformation& Tr, const int nElement,
                                              const mfem::IntegrationPoint& ip,
                                              const int dim) override;
  std::vector<mfem::Coefficient> get_flux_coefficient(mfem::ElementTransformation& Tr,
                                                      const int nElement,
                                                      const mfem::IntegrationPoint& ip) override;

 public:
  ThermalDiffusionFluxNLFormIntegrator(const mfem::ParGridFunction& u_old, const Parameters& params,
                                       std::vector<VARS*> auxvars);
  ~ThermalDiffusionFluxNLFormIntegrator();
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
void ThermalDiffusionFluxNLFormIntegrator<VARS>::get_parameters() {
  DiffusionFluxNLFormIntegrator<VARS>::get_parameters();
}

/**
 * @brief Construct a new
 * ThermalDiffusionFluxNLFormIntegrator<VARS>::ThermalDiffusionFluxNLFormIntegrator object
 *
 * @tparam VARS
 * @param u_old
 * @param params
 * @param auxvars
 */
template <class VARS>
ThermalDiffusionFluxNLFormIntegrator<VARS>::ThermalDiffusionFluxNLFormIntegrator(
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
void ThermalDiffusionFluxNLFormIntegrator<VARS>::check_variables_consistency() {
  std::vector<mfem::ParGridFunction> aux_gf = this->get_aux_gf();
  std::vector<std::vector<std::string>> aux_infos = this->get_aux_infos();

  //==============================================================
  // Get Temperature and thermal diffusivity (aux. variables)
  //==============================================================
  bool temperature_found = false;
  bool diffusivity_found = false;
  for (std::size_t i = 0; i < aux_infos.size(); ++i) {
    const auto& variable_info = aux_infos[i];

    MFEM_VERIFY(!variable_info.empty(), "Empty variable_info encountered.");
    const std::string& symbol = toLowerCase(variable_info.back());

    if (symbol == "Temperature") {
      this->temp_gf_ = aux_gf[i];
      temperature_found = true;
    } else if (symbol == "Diffusivity") {
      this->diffu_gf_ = aux_gf[i];
      diffusivity_found = true;
    }

    if (temperature_found && diffusivity_found) break;
  }

  MFEM_VERIFY(temperature_found,
              "ThermalDiffusionFluxNLFormIntegrator<VARS>::check_variables_consistency: "
              "Temperature not found. Please check your data.");

  MFEM_VERIFY(diffusivity_found,
              "ThermalDiffusionFluxNLFormIntegrator<VARS>::check_variables_consistency: Thermal "
              "diffusivity not found. Please check your data.");
}

/**
 * @brief Destroy the ThermalDiffusionFluxNLFormIntegrator<
 * VARS>::ThermalDiffusionFluxNLFormIntegrator object
 *
 * @tparam VARS
 */
template <class VARS>
ThermalDiffusionFluxNLFormIntegrator<VARS>::~ThermalDiffusionFluxNLFormIntegrator() {}

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
  bool diffusivity_found_{false};
  void check_variables_consistency();

 protected:
  mfem::ParGridFunction temp_gf_;
  mfem::ParGridFunction diffu_gf_;

  void get_parameters() override;
  std::vector<mfem::Vector> get_flux_gradient(mfem::ElementTransformation& Tr, const int nElement,
                                              const mfem::IntegrationPoint& ip,
                                              const int dim) final;

  std::vector<double> get_flux_coefficient(const int nElement,
                                           const mfem::IntegrationPoint& ip) override;

 public:
  ThermalDiffusionFluxNLFormIntegrator(const std::vector<mfem::ParGridFunction>& u_old,
                                       const Parameters& params, std::vector<VARS*> auxvars);
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
    const std::vector<mfem::ParGridFunction>& u_old, const Parameters& params,
    std::vector<VARS*> auxvars)
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

    size_t vsize = variable_info.size();

    MFEM_VERIFY(vsize > 1,
                "ThermalDiffusionFluxNLFormIntegrator<VARS>::check_variables_consistency: at least "
                "two additionnal information are expected for auxiliary variables associated with "
                "this integrator");
    const std::string& symbol = toUpperCase(variable_info.back());

    if (toUpperCase(symbol) == "T") {
      this->temp_gf_ = aux_gf[i];
      temperature_found = true;
    } else if (toLowerCase(symbol) == "lambda") {
      // Diffusivity can be directly supplied within this integrator or overloaded by considering a
      // child class.
      this->diffu_gf_ = aux_gf[i];
      this->diffusivity_found_ = true;
    }

    if (temperature_found && diffusivity_found) break;
  }

  MFEM_VERIFY(temperature_found,
              "ThermalDiffusionFluxNLFormIntegrator<VARS>::check_variables_consistency: "
              "Temperature not found. Please check your data.");
}

/**
 * @brief Return the coefficient part of the thermal diffusion flux
 *
 * @tparam VARS
 * @param Tr
 * @param nElement
 * @param ip
 * @return std::vector<double>
 */
template <class VARS>
std::vector<double> ThermalDiffusionFluxNLFormIntegrator<VARS>::get_flux_coefficient(
    const int nElement, const mfem::IntegrationPoint& ip) {
  std::vector<double> coefficient;
  if (this->diffusivity_found_) {
    coefficient.emplace_back(this->diffu_gf_.GetValue(nElement, ip));
  }
  return coefficient;
}

/**
 * @brief Return the gradient part of the thermal diffusion flux
 *
 * @tparam VARS
 * @param Tr
 * @param nElement
 * @param ip
 * @param dim
 * @return std::vector<mfem::Vector>
 */
template <class VARS>
std::vector<mfem::Vector> ThermalDiffusionFluxNLFormIntegrator<VARS>::get_flux_gradient(
    mfem::ElementTransformation& Tr, const int nElement, const mfem::IntegrationPoint& ip,
    const int dim) {
  std::vector<mfem::Vector> gradient;
  mfem::Vector gradT;
  gradT.SetSize(dim);
  auto temp = SlothGridFunction(this->temp_gf_);
  temp.GetGradient(Tr, this->gradPsi, gradT);

  gradient.emplace_back(gradT);
  return gradient;
}

/**
 * @brief Destroy the ThermalDiffusionFluxNLFormIntegrator<
 * VARS>::ThermalDiffusionFluxNLFormIntegrator object
 *
 * @tparam VARS
 */
template <class VARS>
ThermalDiffusionFluxNLFormIntegrator<VARS>::~ThermalDiffusionFluxNLFormIntegrator() {}

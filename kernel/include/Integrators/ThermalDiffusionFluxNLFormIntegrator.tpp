/**
 * @file ThermalDiffusionFluxNLFormIntegrator.tpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief  Class dedicated to the VF of thermal diffusion flux used in energy balance equation
 * @version 0.1
 * @date 2025-09-05
 *
 * @copyright CEA (C) 2025
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
#pragma once

#include <string>
#include <vector>

#include "Coefficients/AxiCylindricalCoefficient.hpp"
#include "Integrators/SlothGridFunction.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameters.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Retrieve and initialize model parameters from the parameter set.
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 */
template <class VARS>
void ThermalDiffusionFluxNLFormIntegrator<VARS>::get_parameters() {
  DiffusionFluxNLFormIntegrator<VARS>::get_parameters();
}

/**
 * @brief Construct a new ThermalDiffusionFluxNLFormIntegrator object.
 *
 * This constructor initializes the nonlinear form integrator. It forwards the provided previous
 * solution fields, simulation parameters, auxiliary variables, and coefficients to the base
 * DiffusionFlux nonlinear form integrator.
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 *
 * @param u_old        Vector of previous-time-step solution fields.
 * @param params       Paramters that can be used with the integrator.
 * @param auxvars      Auxiliary variables required by the inetgrator.
 * @param coefficients List of coefficients defining material properties.
 *
 */
template <class VARS>
ThermalDiffusionFluxNLFormIntegrator<VARS>::ThermalDiffusionFluxNLFormIntegrator(
    Geometry geometry, const std::vector<mfem::ParGridFunction>& u_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients)
    : DiffusionFluxNLFormIntegrator<VARS>(geometry, u_old, params, auxvars, coefficients) {
  this->check_variables_consistency();
}

/**
 * @brief  Check variables consistency
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
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
 * @brief Compute the diffusion coefficients at the integration point
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param Tr Element transformation.
 * @param nElement number of the element
 * @param ip Integration point
 * @param dim Spatial dimension
 *
 * @return The diffusion coefficients at integration point
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
 * @brief The thermal diffusion flux at integration point
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param Tr Element transformation.
 * @param nElement number of the element
 * @param ip Integration point
 * @param dim Spatial dimension
 *
 * @return The thermal diffusion flux at integration point
 */
template <class VARS>
std::vector<mfem::Vector> ThermalDiffusionFluxNLFormIntegrator<VARS>::get_flux_gradient(
    mfem::ElementTransformation& Tr, [[maybe_unused]] const int nElement,
    [[maybe_unused]] const mfem::IntegrationPoint& ip, const int dim) {
  std::vector<mfem::Vector> gradient;
  mfem::Vector gradT;
  gradT.SetSize(dim);
  auto temp = SlothGridFunction(this->temp_gf_);
  temp.GetGradient(Tr, this->gradPsi, gradT);

  gradient.emplace_back(gradT);
  return gradient;
}

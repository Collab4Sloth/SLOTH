/**
 * @file MassDiffusionFluxNLFormIntegrator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Class dedicated to the VF of  diffusion flux (gradient of chemical potential) used in mass
 * balance equation
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
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "Integrators/SlothGridFunction.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameters.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

/**
 * @brief  Class dedicated to the VF of  diffusion flux (gradient of chemical potential) used in
 * mass balance equation
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 */
template <class VARS>
class MassDiffusionFluxNLFormIntegrator : public DiffusionFluxNLFormIntegrator<VARS> {
 private:
  bool dmu_found_{false};
  bool mu_found_{false};
  bool mobilities_found_{false};
  bool enable_diffusion_chemical_potentials_{false};

  void check_variables_consistency();
  std::vector<mfem::Vector> get_flux_gradient_dmu(mfem::ElementTransformation& Tr,
                                                  const int nElement,
                                                  const mfem::IntegrationPoint& ip, const int dim);
  std::vector<mfem::Vector> get_flux_gradient_mu(mfem::ElementTransformation& Tr,
                                                 const int nElement,
                                                 const mfem::IntegrationPoint& ip, const int dim);
  void get_cross_thermal_flux(mfem::Vector& grad_pot, const mfem::ParGridFunction& potential,
                              mfem::ElementTransformation& Tr, const int nElement,
                              const mfem::IntegrationPoint& ip, const int dim);

 protected:
  std::map<std::string, mfem::ParGridFunction> mu_gf_;
  std::map<std::string, mfem::ParGridFunction> dmu_gf_;
  std::map<std::string, mfem::ParGridFunction> mob_gf_;
  std::vector<mfem::ParGridFunction> temp_gf_;
  bool scale_variables_by_temperature_{false};
  bool scale_coefficients_by_temperature_{false};

  void get_parameters() override;
  std::vector<mfem::Vector> get_flux_gradient(mfem::ElementTransformation& Tr, const int nElement,
                                              const mfem::IntegrationPoint& ip,
                                              const int dim) final;

  std::vector<double> get_flux_coefficient(const int nElement,
                                           const mfem::IntegrationPoint& ip) override;

 public:
  MassDiffusionFluxNLFormIntegrator(const std::vector<mfem::ParGridFunction>& u_old,
                                    const std::vector<mfem::ParGridFunction>& aux_old,
                                    const Parameters& params, std::vector<VARS*> auxvars,
                                    const std::vector<Coefficients>& coefficients);
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Retrieve and initialize model parameters from the parameter set.
 *
 * This function reads the parameters required by the integrator from `params_` and stores them
 * in the corresponding member variables.
 *
 * The following parameters are queried:
 * - `ScaleVariablesByTemperature` (bool, optional): flag to indicate if variables are divided by
 * temperature.
 * - `ScaleCoefficientsByTemperature` (bool, optional): flag to indicate if diffusion coefficients
 * are divided by temperature.
 * - `EnableDiffusionChemicalPotentials` (bool, optional): falg to indicate if diffusion chemical
 * potentials are considered.
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 *
 * @pre `params_` must be initialized and contain the required entries.
 */
template <class VARS>
void MassDiffusionFluxNLFormIntegrator<VARS>::get_parameters() {
  DiffusionFluxNLFormIntegrator<VARS>::get_parameters();

  // Scaling interdiffusion coefficient by RT
  this->scale_variables_by_temperature_ =
      this->params_.template get_param_value_or_default<bool>("ScaleVariablesByTemperature", false);
  this->scale_coefficients_by_temperature_ =
      this->params_.template get_param_value_or_default<bool>("ScaleCoefficientsByTemperature",
                                                              false);
  this->enable_diffusion_chemical_potentials_ =
      this->params_.template get_param_value_or_default<bool>("EnableDiffusionChemicalPotentials",
                                                              false);
}

/**
 * @brief Construct a new MassDiffusionFluxNLFormIntegrator object.
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
MassDiffusionFluxNLFormIntegrator<VARS>::MassDiffusionFluxNLFormIntegrator(
    const std::vector<mfem::ParGridFunction>& u_old,
    const std::vector<mfem::ParGridFunction>& aux_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients)
    : DiffusionFluxNLFormIntegrator<VARS>(u_old, aux_old, params, auxvars, coefficients) {
  this->get_parameters();
  this->check_variables_consistency();
}

/**
 * @brief  Check variables consistency
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 */
template <class VARS>
void MassDiffusionFluxNLFormIntegrator<VARS>::check_variables_consistency() {
  //==========================================================
  // Get chemical potentials and mobilities (aux. variables)
  //==========================================================
  bool temperature_found = false;
  std::set<std::string> set_mob, set_pot, set_dpot;
  for (std::size_t i = 0; i < this->aux_infos_.size(); ++i) {
    const auto& variable_info = this->aux_infos_[i];

    MFEM_VERIFY(!variable_info.empty(), "Empty variable_info encountered.");
    size_t vsize = variable_info.size();

    MFEM_VERIFY(vsize > 0,
                "MassDiffusionFluxNLFormIntegrator<VARS>::check_variables_consistency: at least "
                "one additionnal information are expected for auxiliary variables associated with "
                "this integrator");
    const std::string& symbol = toLowerCase(variable_info.back());

    if (symbol == "mu") {
      MFEM_VERIFY(
          vsize == 2,
          "MassDiffusionFluxNLFormIntegrator<VARS>::check_variables_consistency: error while "
          "getting chemical potentials. Expected [element, 'mu']");

      const std::string& elem_name = variable_info[0];
      this->mu_gf_.emplace(toUpperCase(elem_name), std::move(this->aux_gf_[i]));
      this->mu_found_ = true;
      set_pot.insert(elem_name);
    } else if (symbol == "dmu" && this->enable_diffusion_chemical_potentials_) {
      MFEM_VERIFY(
          vsize == 2,
          "MassDiffusionFluxNLFormIntegrator<VARS>::check_variables_consistency: error while "
          "getting diffusion chemical potentials. Expected [element, 'dmu']");

      const std::string& elem_name = variable_info[0];
      this->dmu_gf_.emplace(toUpperCase(elem_name), std::move(this->aux_gf_[i]));
      this->dmu_found_ = true;
      set_dpot.insert(elem_name);
    } else if (symbol == "inter_mob") {
      // Mobilities can be directly supplied within this integrator or overloaded by considering a
      // child class.
      MFEM_VERIFY(
          vsize >= 2,
          "MassDiffusionFluxNLFormIntegrator<VARS>::check_variables_consistency: error while "
          "getting mobilities.Expected [elem_name, 'mob'].");

      const std::string& elem_name = toUpperCase(variable_info[vsize - 2]);

      this->mob_gf_.emplace(toUpperCase(elem_name), std::move(this->aux_gf_[i]));
      this->mobilities_found_ = true;
      set_mob.insert(elem_name);
    } else if (toUpperCase(variable_info.back()) == "T") {
      // this->temp_gf_ = std::move(aux_gf[i]);
      // const std::string& elem_name = toUpperCase(variable_info[vsize - 2]);
      this->temp_gf_.emplace_back(std::move(this->aux_gf_[i]));
      temperature_found = true;
    }
  }

  MFEM_VERIFY(
      this->dmu_found_ || this->mu_found_,
      "MassDiffusionFluxNLFormIntegrator<VARS>::check_variables_consistency: "
      "Neither chemical potentials, nor diffusion chemical potentials found. At least one of them "
      "is required.");

  if (this->enable_diffusion_chemical_potentials_) {
    MFEM_VERIFY(this->dmu_found_,
                "MassDiffusionFluxNLFormIntegrator<VARS>::check_variables_consistency: "
                "Expected diffusion chemical potentials.");

    if (this->mobilities_found_) {
      MFEM_VERIFY(this->mob_gf_.size() == this->dmu_gf_.size(),
                  "MassDiffusionFluxNLFormIntegrator<VARS>::check_variables_consistency: "
                  "As many mobilities as diffusion chemical potentials are expected.");

      MFEM_VERIFY(set_mob == set_dpot,
                  "MassDiffusionFluxNLFormIntegrator<VARS>::check_variables_consistency: "
                  "Elements must be the same for mobilities and diffusion chemical potentials.");
    }

  } else {
    MFEM_VERIFY(this->mu_found_,
                "MassDiffusionFluxNLFormIntegrator<VARS>::check_variables_consistency: "
                "Expected chemical potentials.");
    if (this->mobilities_found_) {
      MFEM_VERIFY(this->mob_gf_.size() == this->mu_gf_.size(),
                  "MassDiffusionFluxNLFormIntegrator<VARS>::check_variables_consistency: "
                  "As many mobilities as chemical potentials are expected.");

      MFEM_VERIFY(set_mob == set_pot,
                  "MassDiffusionFluxNLFormIntegrator<VARS>::check_variables_consistency: "
                  "Elements must be the same for mobilities and chemical potentials.");
    }
  }
  if (this->scale_coefficients_by_temperature_ || this->scale_variables_by_temperature_) {
    MFEM_VERIFY(
        temperature_found,
        "MassDiffusionFluxNLFormIntegrator<VARS>::check_variables_consistency: "
        "No temperature found. Expected temperature when scaling by temperature is required. ");
  }
}

/**
 * @brief Compute the diffusion flux at the integration point
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param Tr Element transformation.
 * @param nElement number of the element
 * @param ip Integration point
 * @param dim Spatial dimension
 *
 * @return The diffusion flux at integration point
 */
template <class VARS>
std::vector<mfem::Vector> MassDiffusionFluxNLFormIntegrator<VARS>::get_flux_gradient(
    mfem::ElementTransformation& Tr, const int nElement, const mfem::IntegrationPoint& ip,
    const int dim) {
  if (this->enable_diffusion_chemical_potentials_) {
    return this->get_flux_gradient_dmu(Tr, nElement, ip, dim);
  } else {
    return this->get_flux_gradient_mu(Tr, nElement, ip, dim);
  }
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
std::vector<double> MassDiffusionFluxNLFormIntegrator<VARS>::get_flux_coefficient(
    const int nElement, const mfem::IntegrationPoint& ip) {
  std::vector<double> coefficient;
  double scaling_factor = 1.;
  if (this->scale_coefficients_by_temperature_) {
    scaling_factor = Physical::R * this->temp_gf_[0].GetValue(nElement, ip);
  }

  if (this->mobilities_found_) {
    coefficient.reserve(this->mob_gf_.size());
    for (const auto& [component, mob_gf] : this->mob_gf_) {
      coefficient.emplace_back(mob_gf.GetValue(nElement, ip) / scaling_factor);
    }
  }

  return coefficient;
}

/**
 * @brief The diffusion flux calculated with diffusion chemical potentials (dmu)
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param Tr Element transformation.
 * @param nElement number of the element
 * @param ip Integration point
 * @param dim Spatial dimension
 *
 * @return The diffusion flux calculated with diffusion chemical
 * potentials at integration point
 */
template <class VARS>
std::vector<mfem::Vector> MassDiffusionFluxNLFormIntegrator<VARS>::get_flux_gradient_dmu(
    mfem::ElementTransformation& Tr, const int nElement, const mfem::IntegrationPoint& ip,
    const int dim) {
  std::vector<mfem::Vector> gradient;
  gradient.reserve(this->dmu_gf_.size());
  mfem::Vector grad_mu;
  grad_mu.SetSize(dim);
  for (const auto& [component, dmu_gf] : this->dmu_gf_) {
    auto dmu = SlothGridFunction(dmu_gf);
    dmu.GetGradient(Tr, this->gradPsi, grad_mu);
    if (this->scale_variables_by_temperature_) {
      this->get_cross_thermal_flux(grad_mu, dmu, Tr, nElement, ip, dim);
    }
    gradient.emplace_back(grad_mu);
  }
  return gradient;
}

/**
 * @brief The diffusion flux calculated with chemical potentials (mu)
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param Tr Element transformation.
 * @param nElement number of the element
 * @param ip Integration point
 * @param dim Spatial dimension
 *
 * @return The diffusion flux calculated with chemical potentials at integration point
 */
template <class VARS>
std::vector<mfem::Vector> MassDiffusionFluxNLFormIntegrator<VARS>::get_flux_gradient_mu(
    mfem::ElementTransformation& Tr, const int nElement, const mfem::IntegrationPoint& ip,
    const int dim) {
  std::vector<mfem::Vector> gradient;
  gradient.reserve(this->mu_gf_.size());
  mfem::Vector grad_mu;
  grad_mu.SetSize(dim);
  for (const auto& [component, mu_gf] : this->mu_gf_) {
    auto mu = SlothGridFunction(mu_gf);

    mu.GetGradient(Tr, this->gradPsi, grad_mu);
    if (this->scale_variables_by_temperature_) {
      this->get_cross_thermal_flux(grad_mu, mu, Tr, nElement, ip, dim);
    }
    gradient.emplace_back(grad_mu);
  }
  return gradient;
}

/**
 * @brief The thermal diffusion flux at the integration point
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param Tr Element transformation.
 * @param nElement number of the element
 * @param ip Integration point
 * @param dim Spatial dimension
 *
 * @return The thermal diffusion flux at the integration point
 */
template <class VARS>
void MassDiffusionFluxNLFormIntegrator<VARS>::get_cross_thermal_flux(
    mfem::Vector& grad_pot, const mfem::ParGridFunction& potential, mfem::ElementTransformation& Tr,
    const int nElement, const mfem::IntegrationPoint& ip, const int dim) {
  mfem::Vector gradT;
  gradT.SetSize(dim);
  auto temp = SlothGridFunction(this->temp_gf_[0]);
  const auto pot_at_ip = potential.GetValue(nElement, ip);
  const auto temp_at_ip = this->temp_gf_[0].GetValue(nElement, ip);
  const auto dmu_over_square_temp_at_ip = pot_at_ip / (temp_at_ip * temp_at_ip);

  temp.GetGradient(Tr, this->gradPsi, gradT);
  gradT *= dmu_over_square_temp_at_ip;
  grad_pot /= temp_at_ip;
  grad_pot -= gradT;
}

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
  this->get_parameters();

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
  // Get chemical potentials and mobilities (aux. variables)
  //==========================================================
  bool temperature_found = false;
  std::set<std::string> set_mob, set_pot, set_dpot;
  for (std::size_t i = 0; i < aux_infos.size(); ++i) {
    const auto& variable_info = aux_infos[i];

    MFEM_VERIFY(!variable_info.empty(), "Empty variable_info encountered.");
    size_t vsize = variable_info.size();
    const std::string& symbol = toLowerCase(variable_info[vsize - 2]);

    if (symbol == "mu") {
      MFEM_VERIFY(
          variable_info.size() == 3,
          "MassDiffusionFluxNLFormIntegrator<VARS>::check_variables_consistency: error while "
          "getting chemical potentials. Expected [element, 'mu']");

      const std::string& elem_name = variable_info[0];
      this->mu_gf_.emplace(toUpperCase(elem_name), std::move(aux_gf[i]));
      this->mu_found_ = true;
      set_pot.insert(elem_name);
    } else if (symbol == "dmu" && this->enable_diffusion_chemical_potentials_) {
      MFEM_VERIFY(
          variable_info.size() == 3,
          "MassDiffusionFluxNLFormIntegrator<VARS>::check_variables_consistency: error while "
          "getting diffusion chemical potentials. Expected [element, 'dmu']");

      const std::string& elem_name = variable_info[0];
      this->dmu_gf_.emplace(toUpperCase(elem_name), std::move(aux_gf[i]));
      this->dmu_found_ = true;
      set_dpot.insert(elem_name);
    } else if (symbol == "inter_mob") {
      // Mobilities can be directly supplied within this integrator or overloaded by considering a
      // child class.
      MFEM_VERIFY(
          variable_info.size() >= 3,
          "MassDiffusionFluxNLFormIntegrator<VARS>::check_variables_consistency: error while "
          "getting mobilities.Expected [elem_name, 'mob'].");

      const std::string& elem_name = toUpperCase(variable_info[variable_info.size() - 3]);

      this->mob_gf_.emplace(toUpperCase(elem_name), std::move(aux_gf[i]));
      this->mobilities_found_ = true;
      set_mob.insert(elem_name);
    } else if (toUpperCase(variable_info[vsize - 2]) == "T") {
      // this->temp_gf_ = std::move(aux_gf[i]);
      const std::string& elem_name = toUpperCase(variable_info[variable_info.size() - 3]);
      this->temp_gf_.emplace_back(std::move(aux_gf[i]));
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
 * @brief Return the gradient part of the mass diffusion flux
 *
 * @tparam VARS
 * @return std::vector<mfem::Vector>
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
 * @brief Return the coefficient part of the mass diffusion flux
 *
 * @tparam VARS
 * @param Tr
 * @param nElement
 * @param ip
 * @return std::vector<double>
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
 * @brief Return the gradient part of the diffusion flux calculated with diffusion chemical
 * potentials (dmu)
 *
 * @tparam VARS
 * @return std::vector<mfem::Vector>
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
 * @brief Return the gradient part of the diffusion flux calculated with chemical potentials (mu)
 *
 * @tparam VARS
 * @return std::vector<mfem::Vector>
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
 * @brief Return the cross thermal flux
 *
 * @tparam VARS
 * @param potential
 * @return mfem::Vector
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

/**
 * @brief Destroy the MassDiffusionFluxNLFormIntegrator<
 * VARS>::MassDiffusionFluxNLFormIntegrator object
 *
 * @tparam VARS
 */
template <class VARS>
MassDiffusionFluxNLFormIntegrator<VARS>::~MassDiffusionFluxNLFormIntegrator() {}

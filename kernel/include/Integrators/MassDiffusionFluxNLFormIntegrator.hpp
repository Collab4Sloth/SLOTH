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

#pragma once
#include <algorithm>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "Integrators/DiffusionFluxNLFormIntegrator.hpp"
#include "Integrators/SlothGridFunction.hpp"
#include "Options/Options.hpp"
#include "Parameters/Parameters.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

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
                                              const int dim) override final;

  std::vector<double> get_flux_coefficient(const int nElement,
                                           const mfem::IntegrationPoint& ip) override;

 public:
  MassDiffusionFluxNLFormIntegrator(const std::vector<mfem::ParGridFunction>& u_old,
                                    const std::vector<mfem::ParGridFunction>& aux_old,
                                    const Parameters& params, std::vector<VARS*> auxvars,
                                    const std::vector<Coefficients>& coefficients);
};
#include "Integrators/MassDiffusionFluxNLFormIntegrator.tpp"

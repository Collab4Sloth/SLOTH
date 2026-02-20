/**
 * @file MeltingConstantNLFormIntegrator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief NonlinearFormIntegrator for ad-hoc constant melting term
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
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Integrators/MeltingBaseNLFormIntegrator.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Parameters/Parameter.hpp"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

/**
 * @brief Class dedicated to the VF of an ad-hoc constant melting contribution found in
 * Allen-Cahn-type equations.
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 */
template <class VARS>
class MeltingConstantNLFormIntegrator : public MeltingBaseNLFormIntegrator<VARS> {
 private:
  double alpha_;

  void get_parameters();

 protected:
  double get_phase_change_at_ip(mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ir,
                                unsigned int blk, const double u, const double un) override;

 public:
  MeltingConstantNLFormIntegrator(const std::vector<mfem::ParGridFunction>& u_old,
                                  const std::vector<mfem::ParGridFunction>& aux_old,
                                  const Parameters& params, std::vector<VARS*> auxvars,
                                  const std::vector<Coefficients>& coefficients);

  virtual ~MeltingConstantNLFormIntegrator() = default;
};
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Construct a new MeltingConstantNLFormIntegrator object.
 *
 * This constructor initializes the nonlinear form integrator. It forwards the provided previous
 * solution fields, simulation parameters, auxiliary variables, and coefficients to the base melting
 * nonlinear form integrator.
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
MeltingConstantNLFormIntegrator<VARS>::MeltingConstantNLFormIntegrator(
    const std::vector<mfem::ParGridFunction>& u_old,
    const std::vector<mfem::ParGridFunction>& aux_old, const Parameters& params,
    std::vector<VARS*> auxvars, const std::vector<Coefficients>& coefficients)
    : MeltingBaseNLFormIntegrator<VARS>(u_old, aux_old, params, auxvars, coefficients) {
  this->integrator_name_ = "MeltingConstant";

  this->get_parameters();
}

/**
 * @brief Retrieve and initialize model parameters from the parameter set.
 *
 * This function reads the parameters required by the integrator from `params_` and stores them
 * in the corresponding member variables.
 *
 * The following parameters are queried:
 * - `melting_factor` (double, optional):  Scaling factor applied to the melting contribution.
 * Defaults to `1.0` if not provided.
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 *
 * @pre `params_` must be initialized and contain the required entries.
 */
template <class VARS>
void MeltingConstantNLFormIntegrator<VARS>::get_parameters() {
  this->alpha_ = this->params_.template get_param_value_or_default<double>("melting_factor", 1.);
}

/**
 * @brief Compute the value of the phase change at integration point
 *
 * @remark Written for two phase
 *
 * @tparam VARS Template parameter defining the variables used in the integrator.
 *
 * @param Tr Element transformation.
 * @param ir Integration point
 * @param blk   Index of the block.
 * @param u value of the current solution at ip
 * @param un value of the previous solution at ip
 *
 * @return The computed phase change at integration point
 */
template <class VARS>
double MeltingConstantNLFormIntegrator<VARS>::get_phase_change_at_ip(
    [[maybe_unused]] mfem::ElementTransformation& Tr,
    [[maybe_unused]] const mfem::IntegrationPoint& ir, [[maybe_unused]] unsigned int blk,
    [[maybe_unused]] const double u, [[maybe_unused]] const double un) {
  return this->alpha_;
}

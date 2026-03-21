/**
 * @file MeltingTemperatureNLFormIntegrator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief NonlinearFormIntegrator for ad-hoc melting term based on temperature
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

/**
 * @brief Class dedicated to the VF of an ad-hoc temperature melting contribution found in
 * Allen-Cahn-type equations.
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 */
template <class VARS>
class MeltingTemperatureNLFormIntegrator : public MeltingBaseNLFormIntegrator<VARS> {
 private:
  std::vector<mfem::ParGridFunction> temp_gf_;

  double melting_temperature_;
  double melting_enthalpy_;
  void get_parameters();
  void check_variables_consistency();

 protected:
  double get_phase_change_at_ip(mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ir,
                                unsigned int blk, const double u, const double un) override;

 public:
  MeltingTemperatureNLFormIntegrator(const std::vector<mfem::ParGridFunction>& u_old,
                                     const std::vector<mfem::ParGridFunction>& aux_old,
                                     const Parameters& params, std::vector<VARS*> auxvars,
                                     const std::vector<Coefficients>& coefficients);

  virtual ~MeltingTemperatureNLFormIntegrator() = default;
};

#include "Integrators/MeltingTemperatureNLFormIntegrator.tpp"

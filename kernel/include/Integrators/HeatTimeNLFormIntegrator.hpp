/**
 * @file HeatTimeNLFormIntegrator.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief VF of the time derivative
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
#include <string>
#include <utility>
#include <vector>

#include "Integrators/SlothNLFormIntegrator.hpp"
#include "Integrators/TimeNLFormIntegrator.hpp"
#include "MAToolsProfiling/MATimersAPI.hxx"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

/**
 * @brief Class dedicated to the VF of the time-derivativefor heat transfer equation
 *
 * @tparam VARS Template parameter defining the variables used
 *              in the integrator.
 */
template <class VARS>
class HeatTimeNLFormIntegrator : public TimeNLFormIntegrator<VARS> {
 protected:
  void get_coefficients() override;

 public:
  HeatTimeNLFormIntegrator(const std::vector<mfem::ParGridFunction> u_old,
                           const std::vector<mfem::ParGridFunction>& aux_old,
                           const Parameters& params, std::vector<VARS*> auxvars,
                           const std::vector<Coefficients>& coefficients);
  virtual ~HeatTimeNLFormIntegrator() = default;
};

#include "Integrators/HeatTimeNLFormIntegrator.tpp"

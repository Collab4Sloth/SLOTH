/**
 * @file IntegratorsOptions.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Options for Integrators
 * @version 0.1
 * @date 2025-12-12
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

#include <string>

#include "Utils/Utils.hpp"

#pragma once

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

struct Integrators {
  enum value {
    Fick,
    Fourier,
    CahnHilliard,
    AllenCahn,
    TimeDerivative,
    HeatTimeDerivative,
    SplitTimeDerivative,
    MeltingTemperature,MassFlux
  };
  static value from(const std::string&);
};
Integrators::value Integrators::from(const std::string& v) {
  static PhaseFieldPrivate::mmap<Integrators::value> m{
      {"TimeDerivative", Integrators::TimeDerivative},
      {"HeatTimeDerivative", Integrators::HeatTimeDerivative},
      {"SplitTimeDerivative", Integrators::SplitTimeDerivative},
      {"Fick", Integrators::Fick},
      {"Fourier", Integrators::Fourier},
      {"CahnHilliard", Integrators::CahnHilliard},
      {"AllenCahn", Integrators::AllenCahn},
      {"MeltingTemperature", Integrators::MeltingTemperature},
      {"MassFlux", Integrators::MassFlux}};
  return m.find("Integrators", v);
}

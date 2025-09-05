/**
 * @file AnalyticalFunctionOptions.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Options for Analytical functions
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

#include <string>

#include "Utils/Utils.hpp"

#pragma once

struct AnalyticalFunctionsType {
  enum value { Heaviside, Sinusoide, Sinusoide2, HyperbolicTangent, Parabolic, Uniform };
  static value from(const std::string&);
};

AnalyticalFunctionsType::value AnalyticalFunctionsType::from(const std::string& v) {
  static PhaseFieldPrivate::mmap<AnalyticalFunctionsType::value> m{
      {"Heaviside", AnalyticalFunctionsType::Heaviside},
      {"Sinusoide", AnalyticalFunctionsType::Sinusoide},
      {"Sinusoide2", AnalyticalFunctionsType::Sinusoide2},
      {"HyperbolicTangent", AnalyticalFunctionsType::HyperbolicTangent},
      {"Parabolic", AnalyticalFunctionsType::Parabolic},
      {"Uniform", AnalyticalFunctionsType::Uniform}};
  return m.find("AnalyticalFunctionsType", v);
}

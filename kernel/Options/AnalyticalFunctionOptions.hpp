/**
 * @file AnalyticalFunctionOptions.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Options for Analytical functions
 * @version 0.1
 * @date 2025-01-09
 *
 * Copyright CEA (c) 2025
 *
 */

#include <string>

#include "Utils/Utils.hpp"

#pragma once

struct AnalyticalFunctionsType {
  enum value { Heaviside, Sinusoide, Sinusoide2, HyperbolicTangent, Parabolic, Uniform };
  static value from(const std::string&);
};

inline AnalyticalFunctionsType::value AnalyticalFunctionsType::from(const std::string& v) {
  static PhaseFieldPrivate::mmap<AnalyticalFunctionsType::value> m{
      {"Heaviside", AnalyticalFunctionsType::Heaviside},
      {"Sinusoide", AnalyticalFunctionsType::Sinusoide},
      {"Sinusoide2", AnalyticalFunctionsType::Sinusoide2},
      {"HyperbolicTangent", AnalyticalFunctionsType::HyperbolicTangent},
      {"Parabolic", AnalyticalFunctionsType::Parabolic},
      {"Uniform", AnalyticalFunctionsType::Uniform}};
  return m.find("AnalyticalFunctionsType", v);
}

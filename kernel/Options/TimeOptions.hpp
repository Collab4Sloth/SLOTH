/**
 * @file TimeOptions.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Options for Time discretization
 * @version 0.1
 * @date 2025-01-09
 *
 * Copyright CEA (c) 2025
 *
 */

#include <string>

#include "Utils/Utils.hpp"

#pragma once

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////
//////// ODE SOLVER
///////////////////////////////////////////////////
struct TimeScheme {
  enum value { EulerImplicit, EulerExplicit, RungeKutta4 };
  static value from(const std::string&);
};

inline TimeScheme::value TimeScheme::from(const std::string& v) {
  static PhaseFieldPrivate::mmap<TimeScheme::value> m{{"EulerImplicit", TimeScheme::EulerImplicit},
                                                      {"EulerExplicit", TimeScheme::EulerExplicit},
                                                      {"RungeKutta4", TimeScheme::RungeKutta4}};
  return m.find("TimeScheme", v);
}

/**
 * @file PhysicalConvergenceOptions.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Options for PhysicalConvergence objects
 * @version 0.1
 * @date 2025-01-09
 *
 * Copyright CEA (c) 2025
 *
 */

#include <string>

#include "Utils/Utils.hpp"

#pragma once

///////////////////////////////////////////////////
//////// CONVERGENCE
///////////////////////////////////////////////////
struct ConvergenceType {
  enum value { RELATIVE_MAX, ABSOLUTE_MAX };
  static value from(const std::string&);
};

/**
 * @file PhaseFieldOptions.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Options for PhaseField problems
 * @version 0.1
 * @date 2025-01-09
 *
 * Copyright CEA (c) 2025
 *
 */

#include "Utils/Utils.hpp"

#pragma once

enum class PhaseChange { Null, Constant, Calphad };
// struct ThermodynamicsPotentials {
//   enum value { W, F, H, X };
//   static value from(const std::string&);
// };

enum class ThermodynamicsPotentials { W, F, H, X, LOG };
enum class ThermodynamicsPotentialDiscretization { Implicit, Explicit, SemiImplicit };

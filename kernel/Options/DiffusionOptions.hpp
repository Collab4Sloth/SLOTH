/**
 * @file DiffusionOptions.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Option for Mass Diffusion problems
 * @version 0.1
 * @date 2025-01-09
 *
 * Copyright CEA (c) 2025
 *
 */

#include "Utils/Utils.hpp"

#pragma once

///////////////////////////////////////////////////
//////// Diffusion
///////////////////////////////////////////////////
enum class DiffusionCoefficients { Linear };
enum class CoefficientDiscretization {
  Implicit,
  Explicit,
};

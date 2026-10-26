/**
 * @file TimeOptions.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Options for Time discretization
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
#include <string>

#include "Utils/Utils.hpp"

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////
//////// Coefficient
///////////////////////////////////////////////////
/**
 * @brief Time treatment of a coefficient depending on the unknowns.
 *
 * - Implicit: evaluated with the current values of the unknowns, i.e. at the end of the time
 *   step for implicit Euler, and at the end of each stage for a Runge-Kutta scheme
 *   (u_stage = x_i + gamma_i dt k_i, see TransientOperator::ImplicitSolve). Its derivatives are
 *   included in the Jacobian.
 * - Explicit: evaluated with the values given at the beginning of the evaluation, i.e. u_n for
 *   implicit or explicit Euler, but the base vector x_i of the stage (and not u_n) for a
 *   Runge-Kutta scheme (SDIRK...). No contribution to the Jacobian.
 * - SemiImplicit: evaluated with both (current and beginning of the evaluation).
 * - Constant: value independent of the unknowns.
 * - File: value read from a file (FileCoefficient), independent of the unknowns.
 *
 * In all cases, a coefficient depending on time is evaluated at the time of the evaluation:
 * end of the step for implicit Euler, beginning for explicit Euler, time of the stage
 * (t + c_i dt) for a Runge-Kutta scheme (TransientOperator::SetTime).
 */
enum Scheme { Implicit, Explicit, SemiImplicit, Constant, File };

///////////////////////////////////////////////////
//////// ODE SOLVER
///////////////////////////////////////////////////
struct TimeScheme {
  enum value { EulerImplicit, SDIRK23, SDIRK33, EulerExplicit, RungeKutta4 };
  static value from(const std::string&);
};

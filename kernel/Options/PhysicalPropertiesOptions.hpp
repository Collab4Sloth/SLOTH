/**
 * @file PhysicalPropertiesOptions.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Options for physical properties
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

#include "Utils/Utils.hpp"

#pragma once

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * Set of fundamental physical constants used by SLOTH
 */
namespace Physical {
const double R = 8.314462618;     // molar gas constant in J mol-1 K-1 ;
const double NA = 6.02214076e23;  // Avogadro constant  mol-1

}  // namespace Physical

///////////////////////////////////////////////////
//////// Properties
///////////////////////////////////////////////////
enum class Property { Conductivity, Density, HeatCapacity, Diffusion, Mobility };

////////////////////////
//// Conductivity
///////////////////////
enum class Conductivity { Constant, Linear };

////////////////////////
//// Density
///////////////////////
enum class Density { Constant, Linear };

////////////////////////
//// HeatCapacity
///////////////////////
enum class HeatCapacity { Constant, Linear };

////////////////////////
//// Diffusion
///////////////////////
enum class Diffusion { Constant, Linear, Log };

////////////////////////
//// Mobility
///////////////////////
enum class Mobility { Constant, Degenerated };

////////////////////////
//// Lambda
///////////////////////
enum class Lambda { Constant };

////////////////////////
//// Omega
///////////////////////
enum class Omega { Constant };

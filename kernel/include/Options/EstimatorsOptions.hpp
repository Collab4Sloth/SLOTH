/**
 * @file EstimatorsOptions.hpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief  Options for error estimators
 * @version 0.1
 * @date 2026-08-05
 *
 * @copyright CEA (C) 2026
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
#include "Utils/Utils.hpp"

///////////////////////////////////////////////////
//////// ESTIMATORS
///////////////////////////////////////////////////
enum class ErrorEstimatorType { L2_ZIENKIEWICZ_ZHU, KELLY };

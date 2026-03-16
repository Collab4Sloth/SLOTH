/**
 * @file TimeOptions.cpp
 * @author Clément Introïni (clement.introini@cea.fr)
 * @brief Options for Time discretization
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

#include "Options/TimeOptions.hpp"

#include <string>

#include "Utils/Utils.hpp"

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

TimeScheme::value TimeScheme::from(const std::string& v) {
  static PhaseFieldPrivate::mmap<TimeScheme::value> m{{"EulerImplicit", TimeScheme::EulerImplicit},
                                                      {"SDIRK23", TimeScheme::SDIRK23},
                                                      {"SDIRK33", TimeScheme::SDIRK33},
                                                      {"EulerExplicit", TimeScheme::EulerExplicit},
                                                      {"RungeKutta4", TimeScheme::RungeKutta4}};
  return m.find("TimeScheme", v);
}

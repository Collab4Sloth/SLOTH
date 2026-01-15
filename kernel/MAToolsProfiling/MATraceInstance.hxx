/**
 * @file MATraceInstance.hxx
 * @author Raphaël Prat (raphael.prat@cea.fr)
 * @brief
 * @version 0.1
 * @date 2025-12-27
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
#pragma once

#include <MAToolsProfiling/MATraceTypes.hxx>

namespace MATools {
namespace MATrace {
MATrace_point& get_ref_MATrace_point() {
  static MATrace_point instance;
  return instance;
}

MATrace_point& get_MATrace_point() {
  static MATrace_point instance;
  return instance;
}

Trace& get_local_MATrace() {
  static Trace instance;
  return instance;
}
}  // namespace MATrace
}  // namespace MATools

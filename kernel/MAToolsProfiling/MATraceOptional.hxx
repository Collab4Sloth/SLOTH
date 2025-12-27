/**
 * @file MATraceOptional.hxx
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

#include <MAToolsProfiling/MATraceOptional.hxx>
#include <MAToolsProfiling/MAOutput.hxx>
#include <MAToolsProfiling/MATrace.hxx>

namespace MATools {
namespace MATrace {
namespace Optional {
// define default mode values
constexpr bool MATrace_default_mode = false;

extern bool& get_MATrace_mode() {
  static bool _ftm = MATrace_default_mode;
  return _ftm;
}

void active_MATrace_mode() {
  bool& mode = get_MATrace_mode();
  mode = true;
  MATools::MAOutput::printMessage("MATrace_LOG: MATrace is activated");
}

bool is_MATrace_mode() {
  bool ret = get_MATrace_mode();
  return ret;
}
}  // namespace Optional
}  // namespace MATrace
}  // namespace MATools

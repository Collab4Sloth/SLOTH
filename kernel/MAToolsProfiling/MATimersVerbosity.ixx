/**
 * @file MATimersVerbosity.ixx
 * @author Raphaël Prat (raphael.prat@cea.fr)
 * @brief
 * @version 0.1
 * @date 2025-09-08
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
#include <MAToolsProfiling/MAOutput.hxx>
#include <cassert>

namespace MATools {
namespace MATimer {
/**
 * @brief This function displays the name if MATIMERS_VEROBSITY_LEVEL_1 is defined.
 * @param [in] a_name of the chrono section measured.
 * @param [in] a_node_level is the value of the current MATimerNode level.
 */
void print_verbosity_level_1(std::string a_name, int a_node_level) {
#ifdef MATIMERS_VEROBSITY_LEVEL_1
  using namespace MATools::MAOutput;
  assert(a_node_level >= 0);
  assert(a_name != "");
  const std::string level_string = "--";
  std::string message = "Verbosity_message:";
  for (int i = 0; i < a_node_level; i++) {
    message += level_string;
  }
  printMessage(message, ">", a_name);
#endif /* MATIMERS_VEROBSITY_LEVEL_1 */
}
}  // namespace MATimer
}  // namespace MATools

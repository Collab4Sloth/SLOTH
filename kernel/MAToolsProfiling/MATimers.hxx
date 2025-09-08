/**
 * @file MATimers.hxx
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
#include <MAToolsProfiling/EnumTimer.hxx>
#include <MAToolsProfiling/MAOutputManager.hxx>
#include <MAToolsProfiling/MATimerInfo.hxx>
#include <MAToolsProfiling/MATimerNode.hxx>
#include <MAToolsProfiling/Timer.hxx>
#include <cassert>
#include <iostream>

namespace MATools {
namespace MATimer {
bool& is_enable() {
  static bool matimer_enable = false;
  return matimer_enable;
}

void active() {
  auto& enable = is_enable();
  enable = true;
}

void disactive() {
  auto& enable = is_enable();
  enable = false;
}

/**
 * @brief initialize a root timer node. This function has to be followed by finalize function. Do
 * not call this function twice.
 */
void initialize();

/**
 * @brief This function displays each timer node and writes a file with the same information.
 */
void print_and_write_timers();

/**
 * @brief finalize the root timer node. This function has to be called after the intialize function.
 * Do not call this function twice.
 */
void finalize();
}  // namespace MATimer
}  // namespace MATools

#include <MAToolsProfiling/MATimers.ixx>

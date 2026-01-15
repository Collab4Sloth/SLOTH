/**
 * @file MATraceColor.hxx
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

namespace MATools {
namespace MATrace {
struct MATraceRGB{double r; double g; double b; };

constexpr int default_color_size() {
  constexpr int ret = 12;
  return ret;
}

MATraceRGB get_idle_color() {
  MATraceRGB ret = {0.5, 0.5, 0.5};
  return ret;
}

MATraceRGB get_default_color(int color_id) {
  auto size = default_color_size();
  const int i = color_id % size;
  switch (i) {
    case 0 : return {1, 0, 0};
    case 1 : return {0, 0, 1};
    case 2 : return {0, 1, 0};
    case 3 : return {0.5, 0, 0};
    case 4 : return {1, 1, 0};
    case 5 : return {0.5, 0.5, 0};
    case 6 : return {0, 0.5, 0};
    case 7 : return {0, 1, 1};
    case 8 : return {0, 0.5, 0.5};
    case 9 : return {0, 0, 0.5};
    case 10 : return {1, 0, 1};
    case 11 : return {0.5, 0, 0.5};
    default : return {0, 0, 0};
  }
}
}  // namespace MATrace
}  // namespace MATools

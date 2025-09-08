/**
 * @file MADebug.hxx
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

#include "MAToolsProfiling/MATimerNode.hxx"

namespace MATools {

namespace MADebug {

using namespace MATools::MAOutputManager;
using namespace MATools::MATimer;
using namespace MATools::MPI;

/**
 * return print data in a from a MATimerNode
 * @see class MATimerNode
 */
template <enumTimer T>
void debug_print() {
  MATimerNode*& ptr = get_MATimer_node<T>();

  auto my_print = [](MATimerNode* a_ptr, int& shift, double a_runtime) {
    a_ptr->print(shift, a_runtime);
  };

  auto sort_comp = [](MATimerNode* a_ptr, MATimerNode* b_ptr) {
    return a_ptr->get_name() > b_ptr->get_name();
  };

  int shift = 0;
  double d = ptr->get_duration();
  // recursive_sorted_call(my_print, sort_comp, ptr, shift, d);
  recursive_call(my_print, ptr, shift, d);
}

template <enumTimer T>
void debug_write() {
  MATimerNode*& ptr = get_MATimer_node<T>();

  auto my_print = [](MATimerNode* a_ptr, int& shift, double a_runtime) {
    a_ptr->print(shift, a_runtime);
  };

  int shift = 0;
  double d = ptr->get_duration();
  // recursive_sorted_call(my_print, sort_comp, ptr, shift, d);
  recursive_call(my_print, ptr, shift, d);
}

/**
 * return print data in a from a MATimerNode of the master process
 * @see class MATimerNode
 */
template <enumTimer T>
void master_debug_print() {
  MATimerNode*& ptr = get_MATimer_node<T>();

  auto my_print = [](MATimerNode* a_ptr, int& shift, double a_runtime) {
    a_ptr->print_local(shift, a_runtime);
  };
  int shift = 0;
  double d = ptr->get_duration();
  if (is_master()) recursive_call(my_print, ptr, shift, d);
}
};  // namespace MADebug
};  // namespace MATools

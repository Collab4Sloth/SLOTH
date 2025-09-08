/**
 * @file MATimers.ixx
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
#include <MAToolsProfiling/MATimers.hxx>
#include <MAToolsProfiling/MATimersFullTreeMode.hxx>

namespace MATools {
namespace MATimer {
/**
 * @brief initialize a root timer node. This function has to be followed by finalize function. Do
 * not call this function twice.
 */
void initialize() {
  active();
  MATools::MAOutput::printMessage("MATimers_LOG: MATimers initialization ");
  MATimerNode*& root_timer_ptr = MATools::MATimer::get_MATimer_node<ROOT>();
  MATimerNode*& current = MATools::MATimer::get_MATimer_node<CURRENT>();
  assert((root_timer_ptr == nullptr || root_timer_ptr == current) &&
         "MATimer::initialize has already be called");
  root_timer_ptr = new MATimerNode();
  current = root_timer_ptr;
  assert(current != nullptr);
  MATools::MATimer::start_global_timer<ROOT>();
}

/**
 * @brief This function displays each timer node and writes a file with the same information.
 */
void print_and_write_timers() {
  if (is_enable()) {
    MATools::MATimer::end_global_timer<ROOT>();
    MATools::MAOutputManager::write_file();
    MATools::MAOutputManager::print_timetable<ROOT>();
  }
}

/**
 * @brief finalize the root timer node. This function has to be called after the intialize function.
 * Do not call this function twice.
 */
void print() {
  if (is_enable()) {
    using namespace MATools::MATimer::Optional;
    using namespace MATools::MAOutput;
    MATimerNode* root_ptr = MATools::MATimer::get_MATimer_node<ROOT>();
    MATimerNode* current_ptr = MATools::MATimer::get_MATimer_node<CURRENT>();
    assert(root_ptr != nullptr);
    assert(current_ptr != nullptr);

    MATools::MATimer::end_global_timer<ROOT>();

    if (is_print_timetable()) MATools::MAOutputManager::print_timetable<enumTimer::ROOT>();

    if (is_write_file()) {
      MATools::MAOutput::printMessage("MATimers_LOG: Writing timetable ... ");
      MATools::MAOutputManager::write_file();
    }
  }
}

/**
 * @brief finalize the root timer node. This function has to be called after the intialize function.
 * Do not call this function twice.
 */
void finalize() {
  if (is_enable()) {
    using namespace MATools::MATimer::Optional;
    using namespace MATools::MAOutput;
    MATimerNode* root_ptr = MATools::MATimer::get_MATimer_node<ROOT>();
    MATimerNode* current_ptr = MATools::MATimer::get_MATimer_node<CURRENT>();
    root_ptr = nullptr;
    current_ptr = nullptr;
  }
  disactive();
}
};  // namespace MATimer
};  // namespace MATools

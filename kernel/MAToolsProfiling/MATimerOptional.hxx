/**
 * @file MATimerOptional.hxx
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

/**
 * @namespace MATools
 * @brief Namespace containing utility tools for various purposes.
 */
namespace MATools {
/**
 * @namespace MATimer
 * @brief Namespace containing MATimer-related utilities.
 */
namespace MATimer {
/**
 * @namespace Optional
 * @brief Namespace containing optional MATimer configurations.
 */
namespace Optional {
/**
 * @brief Retrieves the reference to the full tree mode configuration.
 * This function retrieves the reference to the full tree mode configuration.
 * @return The reference to the boolean full tree mode.
 */
extern bool& get_full_tree_mode();

/**
 * @brief Retrieves the reference to the print timetable configuration.
 * This function retrieves the reference to the print timetable configuration.
 * @return The reference to the boolean print timetable mode.
 */
extern bool& get_print_timetable();

/**
 * @brief Retrieves the reference to the write file configuration.
 * This function retrieves the reference to the write file configuration.
 * @return The reference to the boolean write file configuration.
 */
extern bool& get_write_file();

/**
 * @brief Activates the full tree mode.
 * This function activates the full tree mode configuration.
 */
void active_full_tree_mode();

/**
 * @brief Enables the print timetable configuration.
 * This function enables the print timetable configuration.
 */
void enable_print_timetable();

/**
 * @brief Enables the write file configuration.
 * This function enables the write file configuration.
 */
void enable_write_file();

/**
 * @brief Disables the print timetable configuration.
 * This function disables the print timetable configuration.
 */
void disable_print_timetable();

/**
 * @brief Disables the write file configuration.
 * This function disables the write file configuration.
 */
void disable_write_file();

/**
 * @brief Checks if the full tree mode is active.
 * This function checks if the full tree mode is active.
 * @return True if the full tree mode is active, false otherwise.
 */
bool is_full_tree_mode();

/**
 * @brief Checks if the print timetable is enabled.
 * This function checks if the print timetable is enabled.
 * @return True if the print timetable is enabled, false otherwise.
 */
bool is_print_timetable();

/**
 * @brief Checks if the write file is enabled.
 * This function checks if the write file is enabled.
 * @return True if the write file is enabled, false otherwise.
 */
bool is_write_file();
}  // namespace Optional
}  // namespace MATimer
}  // namespace MATools

#include <MAToolsProfiling/MATimerOptional.ixx>

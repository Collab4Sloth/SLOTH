#pragma once

/**
 * @namespace MATools
 * @brief Namespace containing utility tools for various purposes.
 */
namespace MATools
{
  /**
   * @namespace MATimer
   * @brief Namespace containing MATimer-related utilities.
   */
  namespace MATimer
  {
    /**
     * @namespace Optional
     * @brief Namespace containing optional MATimer configurations.
     */
    namespace Optional
    {
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
    }
  }
}

#include <MAToolsProfiling/MATimerOptional.ixx>

#include <MAToolsProfiling/MATimerOptional.hxx>
#include <MAToolsProfiling/MAOutput.hxx>

namespace MATools
{
	namespace MATimer
	{
		namespace Optional
		{
			// define some default values
			constexpr bool full_tree_default_mode = false;
			constexpr bool print_timetable_default_mode = true;
			constexpr bool write_file_default_mode = true;

			/*
			 * @brief accessor to the static boolean value of the full tree mode
			 * @return boolean reference of the full tree mode boolean
			 */
			extern bool& get_full_tree_mode()
			{
				static bool _ftm = full_tree_default_mode;
				return _ftm;
			}

			/*
			 * @brief active the full tree mode. Do not use this function if every timers trees have the same hierarchy.
			 */
			void active_full_tree_mode()
			{
				bool& mode = get_full_tree_mode();
				mode = true;
				MATools::MAOutput::printMessage("MATimers_LOG: full tree mode is activated");
			}

			/* 
			 * @brief accessor to check if the full tree mode is activated or not.
			 * @return boolean value. Warning : this boolean is not a reference.
			 */
			bool is_full_tree_mode()
			{
				bool ret = get_full_tree_mode();
				return ret;
			}
			
			/*
			 * @brief accessor to the static boolean value of the print timetable mode
			 * @return boolean reference of the print timetable boolean
			 */
			extern bool& get_print_timetable()
			{
				static bool _pt = print_timetable_default_mode;
				return _pt;
			}

			/*
			 * @brief disable the print timetable mode. 
			 */
			void disable_print_timetable()
			{
				bool& mode = get_print_timetable();
				mode = false;
				MATools::MAOutput::printMessage("MATimers_LOG: print timetable mode is disabled");
			}

			/*
			 * @brief enable the print timetable mode. 
			 */
			void enable_print_timetable()
			{
				bool& mode = get_print_timetable();
				mode = true;
				MATools::MAOutput::printMessage("MATimers_LOG: print timetable mode is activated");
			}

			bool is_print_timetable()
			{
				bool ret = get_print_timetable();
				return ret;
			}
			
			extern bool& get_write_file()
			{
				static bool _wf = write_file_default_mode;
				return _wf;
			}

			/*
			 * @brief disable the print disable mode. 
			 */
			void disable_write_file()
			{
				bool& mode = get_write_file();
				mode = false;
				MATools::MAOutput::printMessage("MATimers_LOG: write file mode is disabled");
			}

			/*
			 * @brief enable the print write file mode. 
			 */
			void enable_write_file()
			{
				bool& mode = get_write_file();
				mode = true;
				MATools::MAOutput::printMessage("MATimers_LOG: write file mode is activated");
			}

			bool is_write_file()
			{
				bool ret = get_write_file();
				return ret;
			}
		}
	}
}

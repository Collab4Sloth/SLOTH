#pragma once
#include <string>

namespace MATools
{
	namespace MATimer
	{
		/**
		 * @brief This function displays the name if MATIMERS_VEROBSITY_LEVEL_1 is defined.
		 * @param [in] a_name of the chrono section measured.
		 * @param [in] a_node_level is the value of the current MATimerNode level.
		 */
		void print_verbosity_level_1(std::string a_name, int a_node_level);
	}
}

#include <MAToolsProfiling/MATimersVerbosity.ixx>

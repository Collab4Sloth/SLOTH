#pragma once

#include <iostream>
#include <cassert>

#include <MAToolsProfiling/EnumTimer.hxx>
#include <MAToolsProfiling/MATimerNode.hxx>
#include <MAToolsProfiling/Timer.hxx>
#include <MAToolsProfiling/MATimerInfo.hxx>
#include <MAToolsProfiling/MAOutputManager.hxx>

namespace MATools
{
	namespace MATimer
	{
    bool& is_enable()
    {
      static bool matimer_enable = false;
      return matimer_enable;
    }

    void active()
    {
      auto& enable = is_enable();
      enable = true;
    }

    void disactive()
    {
      auto& enable = is_enable();
      enable = false;
    }
   
		/**
		 * @brief initialize a root timer node. This function has to be followed by finalize function. Do not call this function twice.  
		 */
		void initialize();

		/**
		 * @brief This function displays each timer node and writes a file with the same information.
		 */
		void print_and_write_timers();

		/**
		 * @brief finalize the root timer node. This function has to be called after the intialize function. Do not call this function twice.  
		 */
		void finalize();
	}
}

#include <MAToolsProfiling/MATimers.ixx>

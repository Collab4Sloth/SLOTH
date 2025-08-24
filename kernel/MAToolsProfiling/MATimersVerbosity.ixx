#include <MAToolsProfiling/MATimersVerbosity.hxx>
#include <MAToolsProfiling/MAOutput.hxx>
#include <cassert>

namespace MATools
{
  namespace MATimer
  {
    /**
     * @brief This function displays the name if MATIMERS_VEROBSITY_LEVEL_1 is defined.
     * @param [in] a_name of the chrono section measured.
     * @param [in] a_node_level is the value of the current MATimerNode level.
     */
    void print_verbosity_level_1(std::string a_name, int a_node_level)
    {
#ifdef MATIMERS_VEROBSITY_LEVEL_1
      using namespace MATools::MAOutput;
      assert( a_node_level>=0 );
      assert( a_name !="" );
      const std::string level_string = "--";
      std::string message = "Verbosity_message:";
      for(int i=0 ; i<a_node_level ; i++)
      {
  message += level_string;
      }
      printMessage(message,">", a_name);
#endif /* MATIMERS_VEROBSITY_LEVEL_1 */
    }
  }
}


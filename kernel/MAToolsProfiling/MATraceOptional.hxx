#pragma once 

#include <MAToolsProfiling/MATraceOptional.hxx>
#include <MAToolsProfiling/MAOutput.hxx>
#include <MAToolsProfiling/MATrace.hxx>

namespace MATools
{
	namespace MATrace
	{
		namespace Optional
		{
			// define default mode values
			constexpr bool MATrace_default_mode = false;

			extern bool& get_MATrace_mode()
			{
				static bool _ftm = MATrace_default_mode;
				return _ftm;
			}

			void active_MATrace_mode()
			{
				bool& mode = get_MATrace_mode();
				mode = true;
				MATools::MAOutput::printMessage("MATrace_LOG: MATrace is activated");
			}

			bool is_MATrace_mode()
			{
				bool ret = get_MATrace_mode();
				return ret;
			}
		};
	};
};

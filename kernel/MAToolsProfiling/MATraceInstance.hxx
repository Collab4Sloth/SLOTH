#pragma once 

#include <MAToolsProfiling/MATraceTypes.hxx>

namespace MATools
{
	namespace MATrace
	{
		MATrace_point& get_ref_MATrace_point()
		{
			static MATrace_point instance;
			return instance;
		}

		MATrace_point& get_MATrace_point()
		{
			static MATrace_point instance;
			return instance;
		}
		
		Trace& get_local_MATrace()
		{
			static Trace instance;
			return instance;
		}
	};
};

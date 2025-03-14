#include "kernel/sloth.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
#include "tests/tests.hpp"

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <onika/log.h>
#include <cassert>

#include <plugin/wrapper.h>

namespace SlothProto
{
	using namespace onika;
	using namespace onika::scg;
		struct SLOTHLog : public OperatorNode
	{
    using SP = Test<3>::SPA;
		ADD_SLOT( Wrapper<SP> , spatial , INPUT ,  DocString{" add doc "});

		inline void execute() override final
		{
      lout << "log" << std::endl;
		}
	};

	// === register factories ===  
	ONIKA_AUTORUN_INIT(log)
	{
		OperatorNodeFactory::instance()->register_factory( "log", make_compatible_operator< SLOTHLog > );
	}
}



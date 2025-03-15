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

	template<typename FECollection, int DIM>
		struct SetTimeDiscretization : public OperatorNode
	{
    using SP = SpatialDiscretization<FECollection,DIM>;
    using VARS = Variables<FECollection, DIM>;
    using PSTCollection = mfem::ParaViewDataCollection;
    using PST = PostProcessing<FECollection, PSTCollection, DIM>;
    using NLFI = AllenCahnNLFormIntegrator<VARS, ThermodynamicsPotentialDiscretization::Implicit,
                                         ThermodynamicsPotentials::W, Mobility::Constant>;
    using OPE = AllenCahnOperator<FECollection, DIM, NLFI>;
    using PB = Problem<OPE, VARS, PST>;
    using CouplingPB = Coupling<PB>;

    ADD_SLOT( double , physical_time , INPUT , REQUIRED , DocString{"physical time of the simulation"});
    ADD_SLOT( double , final_time , INPUT , REQUIRED , DocString{"final time of the simulation"});
    ADD_SLOT( double , dt , INPUT , REQUIRED , DocString{"Time increment"});
    ADD_SLOT( Wrapper<Coupling<PB>>, cc , INPUT , REQUIRED , DocString{" Coupling problem"}); 
    ADD_SLOT( Wrapper<TimeDiscretization<CouplingPB>>, time, OUTPUT , DocString{" Time discretization"}); 

		inline void execute() override final
		{
      auto time_params = Parameters(Parameter("initial_time", *physical_time), Parameter("final_time", *final_time),  Parameter("time_step", *dt));
      time->alias(
        new TimeDiscretization(time_params, cc->get())
      );
		}
	};

	// === register factories ===  
	ONIKA_AUTORUN_INIT(set_time_discretization)
	{
		using setTimeDiscr3D = SetTimeDiscretization<mfem::H1_FECollection,3>;
		OperatorNodeFactory::instance()->register_factory( "set_time_discretization_3d", make_compatible_operator< setTimeDiscr3D > );
	}
}



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
		struct SlothSolve : public OperatorNode
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

    ADD_SLOT( Wrapper<TimeDiscretization<CouplingPB>>, time, INPUT_OUTPUT , DocString{" Time discretization"}); 
		inline void execute() override final
		{
      auto t = time->get();
      t.solve();
		}
	};

	// === register factories ===  
	ONIKA_AUTORUN_INIT(solve)
	{
		using Solve3D = SlothSolve<mfem::H1_FECollection,3>;
		OperatorNodeFactory::instance()->register_factory( "solve_3d", make_compatible_operator< Solve3D > );
	}
}



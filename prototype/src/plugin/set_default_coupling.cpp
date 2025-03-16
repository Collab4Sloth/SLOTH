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
		struct SetAllenCahnDefaultCoupling : public OperatorNode
	{
    using SP = SpatialDiscretization<FECollection,DIM>;
    using VARS = Variables<FECollection, DIM>;
    using PSTCollection = mfem::ParaViewDataCollection;
    using PST = PostProcessing<FECollection, PSTCollection, DIM>;
    using NLFI = AllenCahnNLFormIntegrator<VARS, ThermodynamicsPotentialDiscretization::Implicit,
                                         ThermodynamicsPotentials::W, Mobility::Constant>;
    using OPE = AllenCahnOperator<FECollection, DIM, NLFI>;
    using PB = Problem<OPE, VARS, PST>;

    ADD_SLOT( Wrapper<PB>      , pb , INPUT  , DocString{" Allen Cahn problem"}); 
    ADD_SLOT( Wrapper<Coupling<PB>>, cc , OUTPUT , DocString{" Coupling problem"}); 

		inline void execute() override final
		{
      auto& ACProblem = pb->get();
      cc->wrap(new Coupling<PB>("Default Coupling", ACProblem));
		}
	};

	// === register factories ===  
	ONIKA_AUTORUN_INIT(allen_cahn_default_coupling)
	{
		using setH1ACCC3D = SetAllenCahnDefaultCoupling<mfem::H1_FECollection,3>;
		OperatorNodeFactory::instance()->register_factory( "allen_cahn_h1_default_coupling", make_compatible_operator< setH1ACCC3D > );
	}
}



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
		struct SetAllenCahnProblem : public OperatorNode
	{
    using SP = SpatialDiscretization<FECollection,DIM>;
    using VARS = Variables<FECollection, DIM>;
    using PSTCollection = mfem::ParaViewDataCollection;
    using PST = PostProcessing<FECollection, PSTCollection, DIM>;
    using NLFI = AllenCahnNLFormIntegrator<VARS, ThermodynamicsPotentialDiscretization::Implicit,
                                         ThermodynamicsPotentials::W, Mobility::Constant>;
    using OPE = AllenCahnOperator<FECollection, DIM, NLFI>;
    using PB = Problem<OPE, VARS, PST>;

		ADD_SLOT( double, crit_cvg , INPUT , REQUIRED , DocString{"Convergence criterion"});
		ADD_SLOT( double, mobility , INPUT , REQUIRED , DocString{"Two-phase mobility"});
		ADD_SLOT( double, epsilon  , INPUT , REQUIRED , DocString{"Interface thickness"});
    ADD_SLOT( double, sigma    , INPUT , REQUIRED , DocString{"Interfacial energy"}); 

    ADD_SLOT( std::string   , name    , INPUT , "SlothProblemDefaultName" );
		ADD_SLOT( Wrapper<SP>   , spatial , INPUT , DocString{" add doc "});
    ADD_SLOT( Wrapper<VARS> , vars    , INPUT , DocString{" List of Variable(s)"}); 
    ADD_SLOT( Wrapper<PB>   , pb      , OUTPUT , DocString{" Allen Cahn problem"}); 
    ADD_SLOT( Wrapper<Parameters> , pp_params , INPUT , DocString{"List of parameters used to define the post processings"} );

    // WARNING I need to keep this structure alive, design issue here
    ADD_SLOT( Wrapper<PST>, pst , INPUT_OUTPUT , DocString{"Design issue"}); 
    ADD_SLOT( Wrapper<Parameters>, test_pp , INPUT_OUTPUT , DocString{"Design issue"}); 
		inline void execute() override final
		{
      double mob = *mobility;
      double eps = *epsilon;
      double sig = *sigma;
      auto* spa = spatial->get_ptr();
      auto& ACProblem = *pb;

      // compute some quantities
      const double lambda = 3. * sig * eps / 2.;
      const double omega = 12. * sig / eps;
 
      test_pp->alias(
        new Parameters(
          Parameter("epsilon", eps), 
          Parameter("sigma", sig),
          Parameter("lambda", lambda), 
          Parameter("omega", omega)
        )
      );

      OPE oper(spa, test_pp->get(), TimeScheme::EulerImplicit);
      oper.overload_mobility(Parameters(Parameter("mob", mob))); 
      PhysicalConvergence convergence(ConvergenceType::ABSOLUTE_MAX, *crit_cvg);

      pst->alias(
        new PST(spa, pp_params->get())
      );

      ACProblem.alias(
        new PB(*name, oper, vars->get(), pst->get(), convergence)
      ); 
		}
	};

	// === register factories ===  
	ONIKA_AUTORUN_INIT(generate_hypercube)
	{
		using setH1ACPB3D = SetAllenCahnProblem<mfem::H1_FECollection,3>;
		OperatorNodeFactory::instance()->register_factory( "set_h1_allen_cahn_problem", make_compatible_operator< setH1ACPB3D > );
	}
}



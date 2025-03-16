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
		struct GenerateHypercube : public OperatorNode
	{
    using SP = SpatialDiscretization<FECollection,DIM>;
		ADD_SLOT( std::vector<double>, discretization, INPUT, REQUIRED, DocString{"discretization of the hypercube"});
		ADD_SLOT( std::vector<int>, grid_dim, INPUT, REQUIRED, DocString{"grid dimension of the hypercube"});
		ADD_SLOT( int , refinement_level , INPUT , 0 /*default */, DocString{"Level of Uniforme Refinenemt, default is 0"});
		ADD_SLOT( int , fe_order , INPUT , 1 /*default */, DocString{"Finite Element order, default is 1"});
		ADD_SLOT( Wrapper<SP> , spatial , OUTPUT ,  DocString{" add doc "});
    ADD_SLOT( int , fake , INPUT); 

		inline void execute() override final
		{
			static_assert(DIM == 2 || DIM == 3);
			auto& disc = *discretization;
			auto& gd = *grid_dim;
			assert(DIM == disc.size());
			assert(DIM == gd.size());

			const std::string spatial_name = "InlineSquareWithTetraedres";
			if constexpr( DIM == 2)
			{
				auto description = std::make_tuple(gd[0], gd[1], disc[0], disc[1]);
        spatial->wrap(
          new SP(spatial_name, *fe_order, *refinement_level, description)
        );
			}
			if constexpr( DIM == 3)
			{
				auto description = std::make_tuple(gd[0], gd[1], gd[2], disc[0], disc[1], disc[2]);
        spatial->wrap(
          new SP(spatial_name, *fe_order, *refinement_level, description)
        );
			}
		}
	};

	// === register factories ===  
	ONIKA_AUTORUN_INIT(generate_hypercube)
	{
		using GenerateH1Square = GenerateHypercube<mfem::H1_FECollection,2>;
		using GenerateH1Cube   = GenerateHypercube<mfem::H1_FECollection,3>;
		OperatorNodeFactory::instance()->register_factory( "generate_h1_square", make_compatible_operator< GenerateH1Square > );
		OperatorNodeFactory::instance()->register_factory( "generate_h1_cube", make_compatible_operator< GenerateH1Cube > );
	}
}



#include "kernel/sloth.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
#include "tests/tests.hpp"

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <onika/math/basic_types.h>
#include <onika/log.h>
#include <cassert>

#include <plugin/wrapper.h>

namespace SlothProto
{

	using namespace onika;
	using namespace onika::scg;
	using namespace onika::math;
  struct SetVariablesTest1 : public OperatorNode
	{
    using FECollection = mfem::H1_FECollection;
    static constexpr int DIM = 3;
    using SP = SpatialDiscretization<FECollection,DIM>;
    using BCS = BoundaryConditions<FECollection,DIM>;
    using VARS = Variables<FECollection, DIM>;
		ADD_SLOT( Vec3d, center, INPUT, REQUIRED, DocString{"Center of a sphere"});
		ADD_SLOT( Vec3d, no_idea, INPUT, REQUIRED, DocString{"Looks important"});
		ADD_SLOT( double, thickness, INPUT, REQUIRED, DocString{""});
		ADD_SLOT( double, radius, INPUT, REQUIRED, DocString{"Radius of the sphere"});
    ADD_SLOT( double, epsilon, INPUT, REQUIRED, DocString{"Maybe the solution"});
    ADD_SLOT( std::string, field_name, INPUT, "phi", DocString{" Field name, default is phi"});

		ADD_SLOT( Wrapper<SP> , spatial , INPUT ,  DocString{" add doc "});
    ADD_SLOT( Wrapper<BCS> , bcs , INPUT , DocString{" List of boundary conditions "}); 
    ADD_SLOT( Wrapper<VARS> , vars , OUTPUT , DocString{" List of Variable(s)"}); 

		inline void execute() override final
		{
      static_assert(DIM == 3);
      assert(*radius > 0.);
      assert(*thickness > 0.);

      Vec3d c = *center;
      Vec3d a = *no_idea;

      auto initial_condition = AnalyticalFunctions<DIM>(AnalyticalFunctionsType::HyperbolicTangent, c.x, c.y, c.z, a.x, a.y, a.z, *thickness, *radius);

      auto analytical_solution = AnalyticalFunctions<DIM>(AnalyticalFunctionsType::HyperbolicTangent, c.x, c.y, c.z, a.x, a.y, a.z, *epsilon, *radius);

      auto var = Variable<FECollection,DIM>(spatial->get_ptr(), bcs->get(),*field_name, 2, initial_condition, analytical_solution);
      vars->wrap(
        new VARS(var)
      );
		}
	};

	// === register factories ===  
	ONIKA_AUTORUN_INIT(generate_hypercube)
	{
		OperatorNodeFactory::instance()->register_factory( "set_variable_test1", make_compatible_operator< SetVariablesTest1 > );
	}
}



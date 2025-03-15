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
		struct SetBoundaryConditions : public OperatorNode
	{
    using SP = SpatialDiscretization<FECollection,DIM>;
    using BCS = BoundaryConditions<FECollection,DIM>;
		ADD_SLOT( std::vector<double>, values, INPUT, REQUIRED, DocString{""});
		ADD_SLOT( std::vector<std::string>, types, INPUT, REQUIRED, DocString{""});
		ADD_SLOT( std::vector<std::string>, names, INPUT, REQUIRED, DocString{""});

		ADD_SLOT( Wrapper<SP> , spatial , INPUT ,  DocString{" add doc "});
    ADD_SLOT( Wrapper<BCS> , bcs , OUTPUT , DocString{" add doc "}); 

		inline void execute() override final
		{
      auto& _types = *types;
      auto& _values = *values;
      auto& _names = *names;
      auto* spa = spatial->get_ptr();
      if(_types.size() != _names.size() || _types.size() != _values.size())
      {
        lout << "Error, types, names, and values have to be the same size" << std::endl;
        lout << "types size  = " << _types.size() << std::endl;
        lout << "names size  = " << _names.size() << std::endl;
        lout << "vamues size = " << _values.size() << std::endl;
        std::abort();
      }
      std::vector<Boundary> boundaries;
      for(size_t i = 0 ; i < _names.size() ; i++)
      {
        boundaries.push_back(Boundary(_names[i], i, _types[i], _values[i]));
      }
      auto& boundary_conditions = bcs->get();
      boundary_conditions.Initialize(spa, boundaries); 
		}
	};

	// === register factories ===  
	ONIKA_AUTORUN_INIT(generate_hypercube)
	{
		using setH1BCS = SetBoundaryConditions<mfem::H1_FECollection,3>;
		OperatorNodeFactory::instance()->register_factory( "set_h1_bcs", make_compatible_operator< setH1BCS > );
	}
}



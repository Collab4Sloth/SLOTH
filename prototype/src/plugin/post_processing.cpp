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
		struct SLOTHPostProcessing : public OperatorNode
	{
    ADD_SLOT( int , paraview_frequency , INPUT , REQUIRED );
    ADD_SLOT( int , verbosity , INPUT , REQUIRED );
    ADD_SLOT( std::string , calculation_path , INPUT , "Problem" );
    ADD_SLOT( std::string , main_folder_path , INPUT , "Saves" );
    ADD_SLOT( Wrapper<Parameters> , pp_params , OUTPUT , DocString{"List of parameters used to define the post processings"} );

		inline void execute() override final
		{
      int pf = *paraview_frequency;
      int v  = *verbosity;
      auto& pp = *pp_params;
      lout << " Setting post processing ... " << std::endl;
      lout << " Folder: "    << *main_folder_path << "/" << *calculation_path << std::endl;
      lout << " Frequency: " << pf << std::endl;
      lout << " Level of verbosity: " << v << std::endl;
      pp.wrap(
        new Parameters(
          Parameter("main_folder_path", *main_folder_path), 
          Parameter("calculation_path", *calculation_path), 
          Parameter("frequency", pf), 
          Parameter("level_of_detail", v))
      );
		}
	};

	// === register factories ===  
	ONIKA_AUTORUN_INIT(log)
	{
		OperatorNodeFactory::instance()->register_factory( "post_processing", make_compatible_operator< SLOTHPostProcessing > );
	}
}



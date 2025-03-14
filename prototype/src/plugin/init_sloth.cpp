#include "kernel/sloth.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]
#include "tests/tests.hpp"

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <onika/log.h>

namespace SlothProto
{
  using namespace onika;
  using namespace onika::scg;
  struct InitSloth : public OperatorNode
  {
    ADD_SLOT( bool    , profiling , INPUT_OUTPUT , true , DocString{"Acitivate the profiling"});
    ADD_SLOT( int    , fake ,  OUTPUT);

    inline void execute() override final
    {
//      mfem::Mpi::Init(argc, argv);
      mfem::Hypre::Init();
      if( *profiling )
      {
        Profiling::getInstance().enable();
      }
    }
  };

  // === register factories ===  
  ONIKA_AUTORUN_INIT(init_sloth)
  {
    OperatorNodeFactory::instance()->register_factory( "init_sloth", make_compatible_operator< InitSloth > );
  }
}



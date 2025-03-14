#include <onika/app/api.h>

int main(int argc,char*argv[])
{
  // ============= run simulation ================
  auto app = onika::app::init(argc,argv);
  if( app->get_error_code() >= 0 ) return app->get_error_code();

  // run full simulation graph
  onika::app::run( app );

  // finalize simulation graph
  onika::app::end( app );

  return 0;
}


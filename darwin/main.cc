#include "ce.h"
#include "individual.h"
#include "darwin.h"
#include "print_xmgrace.h"

#include <eo/eoScalarFitness.h>
#include <mpi/mpi_object.h>

int main(int argc, char *argv[]) 
{
  mpi::main( argc, argv );
  darwin::printxmg.init("convex_hull.agr");
  try
  {
    typedef darwin::Individual< CE::Object, eoMinimizingFitness > t_Individual;
    darwin::Darwin< t_Individual, CE::Evaluator > ga;
    if ( not  ga.Load("input.xml") )
      throw ""; 
    ga.run();
  }
  catch ( std::exception &e )
  {
    std::cerr << "Caught error while running lada" << std::endl
	      << e.what();
  }
  return 0;
}

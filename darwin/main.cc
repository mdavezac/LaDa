#ifdef _CE
  #include "ce.h"
  namespace Functional = CE;
#endif
#ifdef _PESCAN
  #include "pescan.h"
  namespace Functional = BandGap;
#endif
#include "individual.h"
#include "darwin.h"
#include "print_xmgrace.h"

#include <eo/eoScalarFitness.h>
#ifdef _MPI
#  include <mpi/mpi_object.h>
#endif

int main(int argc, char *argv[]) 
{
#ifdef _MPI
  mpi::main(argc, argv);
#endif
  darwin::printxmg.init("convex_hull.agr");
  try
  {
    typedef darwin::Individual< Functional::Object, eoMinimizingFitness > t_Individual;
    darwin::Darwin< t_Individual, Functional::Evaluator > ga;
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

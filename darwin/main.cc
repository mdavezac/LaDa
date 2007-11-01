//
//  Version: $Id$
//
#include <stdexcept>       // std::runtime_error
#include <eo/eoScalarFitness.h>

#ifdef _PESCAN
#ifdef _CE 
#undef _CE
#endif
  #include "bandgap.h"
  typedef BandGap :: Evaluator t_Evaluator;
#elif _CE
  #include "groundstate.h"
  typedef GroundState :: Evaluator t_Evaluator;
#elif _MOLECULARITY
  #include "molecularity.h"
  typedef Molecularity :: Evaluator t_Evaluator;
#else 
#error Need to define _CE or _PESCAN or _MOLECULARITY
#endif
#include "individual.h"
#include "darwin.h"
#include "print/xmg.h"

#ifdef _MPI
#  include "mpi/mpi_object.h"
#endif

int main(int argc, char *argv[]) 
{
#ifdef _MPI
  mpi::main(argc, argv);
#endif
  try
  {
    GA::Darwin< t_Evaluator > ga;
    if ( not  ga.Load("input.xml") )
      throw std::runtime_error( "Could not load input!!\n" ); 
    ga.run();
  }
  catch ( std::exception &e )
  {
    std::cerr << "Caught error while running lada" << std::endl
	      << e.what();
  }
  return 0;
}

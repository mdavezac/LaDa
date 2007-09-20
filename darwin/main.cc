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
  typedef BandGap :: Object t_Object;
#elif _CE
  #include "ce.h"
  typedef CE :: Evaluator t_Evaluator;
  typedef CE :: Object t_Object;
#elif _MOLECULARITY
  #include "molecularity.h"
  typedef Molecularity :: Evaluator t_Evaluator;
  typedef Molecularity :: Object t_Object;
#else 
  Need to define _CE or _PESCAN
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
  Print::xmg.init("convex_hull.agr");
  try
  {
    typedef Traits::GA< t_Evaluator > t_GATraits;
    GA::Darwin< t_GATraits > ga;
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

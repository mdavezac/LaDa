//
//  Version: $Id$
//
#include <stdexcept>       // std::runtime_error

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
#elif _EMASS
  #include "emass.h"
  typedef eMassSL :: Evaluator t_Evaluator;
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
  std::string filename("input.xml");
  if( argc > 1 )
  {
    std::ostringstream sstr;
    for( types::t_int i = 1; i < argc; ++i )
      sstr << argv[i] << " "; 
    std::istringstream istr( sstr.str() );
    while ( istr.good() )
    {
      std::string is_op;
      istr >> is_op; is_op = Print::StripEdges( is_op );
      if( is_op.empty() ) continue;
      else if(     istr.good()
               and (is_op == "-i" or is_op == "--input") ) istr >> filename;
      else if( is_op == "-h" or is_op == "--help" )
        std::cout << "Command-line options:\n\t -h, --help this message"
                  << "\n\t -i, --input XML input file (default: input.xml)\n\n";
    }
    filename = Print::reformat_home( filename );
    if( filename != "input.xml" )
      std::cout << "Reading from input file " << filename << std::endl;
  }

  try
  {
    GA::Darwin< t_Evaluator > ga;
    if ( not  ga.Load(filename) )
    {
      std::ostringstream sstr;
      sstr << __FILE__ << ", line: " << __LINE__ << "\n"
           << "Could not load input from file " << filename << "\n";
      throw std::runtime_error( sstr.str() ); 
    }
    ga.run();
  }
  catch ( std::exception &e )
  {
    std::cerr << "Caught error while running lada" << std::endl
	      << e.what();
  }
  return 0;
}

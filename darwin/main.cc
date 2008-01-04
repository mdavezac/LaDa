//
//  Version: $Id$
//
#include <stdexcept>       // std::runtime_error

#ifdef _PESCAN
  #include "bandgap.h"
  typedef BandGap :: Evaluator t_Evaluator;
#define __PROGNAME__ "bandgap_opt"
#elif _CE
  #include "groundstate.h"
  typedef GroundState :: Evaluator t_Evaluator;
#define __PROGNAME__ "ce_opt"
#elif _MOLECULARITY
  #include "molecularity.h"
  typedef Molecularity :: Evaluator t_Evaluator;
#define __PROGNAME__ "layered_opt"
#elif _EMASS
  #include "emass.h"
  typedef eMassSL :: Evaluator t_Evaluator;
#define __PROGNAME__ "emass_opt"
#else 
#error Need to define _CE or _PESCAN or _MOLECULARITY
#endif
#include "individual.h"
#include "darwin.h"
#include "print/xmg.h"

#ifdef _MPI
#  include "mpi/mpi_object.h"
#endif

void parse_cli( int argc, char *argv[], std::string &_filename )
{
  if( argc < 2 ) return;
  
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
             and (is_op == "-i" or is_op == "--input") ) istr >> _filename;
    else if( is_op == "-h" or is_op == "--help" )
      std::cout << "Command-line options:\n\t -h, --help this message"
                << "\n\t -i, --input XML input file (default: input.xml)\n\n";
  }
  _filename = Print::reformat_home( _filename );
}

int main(int argc, char *argv[]) 
{
  std::string filename("input.xml");
#ifdef _MPI
  mpi::main(argc, argv);
  if( mpi::main.is_root_node() )
#endif
  parse_cli( argc, argv, filename );

#ifdef _MPI
  mpi::BroadCast bc(mpi::main);
  if( not mpi::main.is_root_node() ) filename.clear();
  bc << filename 
     << mpi::BroadCast::allocate 
     << filename 
     << mpi::BroadCast::broadcast 
     << filename
     << mpi::BroadCast::clear;
  if( mpi::main.is_root_node() )
#endif
  if( filename != "input.xml" )
    std::cout << "Reading from input file " << filename << std::endl;

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
    std::ostringstream sstr;
#ifdef _MPI
    sstr.str("");
    sstr << "\nProcessor " << mpi::main.rank() + 1
         << " of " << mpi::main.size()
         << " says:\n";
#endif
    sstr << "Caught error while running " << __PROGNAME__ 
         << "\n" << e.what() << "\n";
#ifdef _MPI
    mpi::AllGather bc( mpi::main );
    std::string message = sstr.str();
    bc << mpi::BroadCast::clear
       << message
       << mpi::BroadCast::allocate
       << message
       << mpi::BroadCast::broadcast;
    sstr.str("");
    for( types::t_int n = mpi::main.size(); n > 0; --n )
    {
      bc.serialize( message );
      sstr << message;
    }
    bc << mpi::BroadCast::clear;
    if( mpi::main.is_root_node() ) 
#endif
      std::cerr << sstr.str() << std::endl;

  }
  return 0;
}

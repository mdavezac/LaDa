//
//  Version: $Id$
//
#include <stdexcept>       // std::runtime_error

#include <revision.h>
#ifdef _PESCAN
  #include "bandgap.h"
  typedef BandGap :: Evaluator t_Evaluator;
#define __PROGNAME__ "Band-Gap Optimization"
#elif defined(_CE)
  #include "groundstate.h"
  typedef GroundState :: Evaluator t_Evaluator;
#define __PROGNAME__ "Cluster Expansion Optimization"
#elif defined(_MOLECULARITY)
  #include "molecularity.h"
  typedef Molecularity :: Evaluator t_Evaluator;
#define __PROGNAME__ "Band-Gap Optimization for Epitaxial Structure"
#elif defined(_EMASS)
  #include "emass.h"
  typedef eMassSL :: Evaluator t_Evaluator;
#define __PROGNAME__ "emass_opt"
#else 
#error Need to define _CE or _PESCAN or _MOLECULARITY
#endif
#include "individual.h"
#include "darwin.h"
#include "print/xmg.h"

#include <opt/debug.h>
#include <mpi/mpi_object.h>

bool parse_cli( int argc, char *argv[], std::string &_filename )
{
  bool dostop = false;

  __NOTMPIROOT( goto syncfilename; )

  { // serial and root node 
    std::ostringstream sstr;
    for( types::t_int i = 1; i < argc; ++i )
      sstr << argv[i] << " "; 
    std::istringstream istr( sstr.str() );
    
    while ( istr.good() and (not dostop) )
    {
      std::string is_op;
      istr >> is_op; is_op = Print::StripEdges( is_op );
      if( is_op.empty() ) continue;
      else if(     istr.good()
               and (is_op == "-i" or is_op == "--input") ) istr >> _filename;
      else if( is_op == "-h" or is_op == "--help" )
      {
        std::cout << "\n" << __PROGNAME__ << " from the " << PACKAGE_STRING << " package."
                  << "\nCommand-line options:\n\t -h, --help this message"
                  << "\n\t -v, --version Subversion Revision and Package version"
                  << "\n\t -i, --input XML input file (default: input.xml)\n\n";
        dostop = true;
      }
      else if( is_op == "-v" or is_op == "--version" )
      {
        std::cout << "\n" << __PROGNAME__ << " from the " << PACKAGE_STRING << " package\n"
                  << "Subversion Revision: " << SVN::Revision << "\n\n"; 
        dostop = true;
      }
    }
    _filename = Print::reformat_home( _filename );
  }

syncfilename:
// __TRYMPICODE(
    mpi::BroadCast bc(mpi::main);
    __NOTMPIROOT(  _filename.clear(); )
    bc << _filename  << dostop
       << mpi::BroadCast::allocate 
       << _filename  << dostop
       << mpi::BroadCast::broadcast 
       << _filename  << dostop
       << mpi::BroadCast::clear; // ,
 //  "Error while syncing input filename " << _filename << " between procs.\n"
 //)

dostopnow:
  if( dostop ) return false;
  
  __NOTMPIROOT( return true; )

  if( _filename != "input.xml" )
    std::cout << "Reading from input file " << _filename << std::endl;
  return true;
}

int main(int argc, char *argv[]) 
{
  std::string filename("input.xml");
  __MPICODE( mpi::main(argc, argv); )

  if( not parse_cli( argc, argv, filename ) ) return true;


  try
  {
    GA::Darwin< t_Evaluator > ga;
    __DOASSERT( not ga.Load(filename),
                "Could not load input from file " << filename << "\n" )

    ga.run();
  }
  catch ( std::exception &e )
  {
    std::ostringstream sstr;
    __MPICODE( sstr << "\nProcessor " << mpi::main.rank() + 1
                    << " of " << mpi::main.size()
                    << " says:\n"; )

    sstr << "Caught error while running " << __PROGNAME__ 
         << "\n" << e.what() << "\n";

    __MPICODE( 
      mpi::AllGather bc( mpi::main );
      std::string message = sstr.str();
      bc << mpi::BroadCast::clear
         << message
         << mpi::BroadCast::allocate
         << message
         << mpi::BroadCast::broadcast;
      sstr.str("");
      while( bc.serialize( message ) ) sstr << message;
      bc << mpi::BroadCast::clear;
    )
    __NOTMPIROOT( return 0; )
    std::cerr << sstr.str() << std::endl;
  }
  return 0;
}

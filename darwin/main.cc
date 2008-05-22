//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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

#include <stdexcept>       // std::runtime_error
#include <boost/program_options.hpp>

#include <revision.h>
#include <print/xmg.h>
#include <print/stdout.h>
#include <opt/debug.h>
#include <mpi/mpi_object.h>

#include "individual.h"
#include "darwin.h"

int main(int argc, char *argv[]) 
{
  try
  {
    __MPICODE( boost::mpi::environment env(argc, argv); )
    namespace po = boost::program_options;
 
    std::string filename("input.xml");
    __MPICODE( mpi::main(argc, argv); )
    unsigned reruns = 1;
 
    po::options_description generic("Generic Options");
    generic.add_options()
           ("help,h", "produces this help message.")
           ("version,v", "prints version string.");
    po::options_description specific("GA Options");
    specific.add_options()
        ("input,i", po::value<std::string>()->default_value("input.xml"), 
         "input filename." )
        ("reruns,r", po::value<unsigned>()->default_value(1),
                     "number of times to run the algorithm.\n"
                     "Is equivalent to manually re-launching the program.\n");
 
    po::options_description all;
    all.add(generic).add(specific);
 
    po::positional_options_description p;
    p.add("input", 1);
 
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(all).positional(p).run(), vm);
    po::notify(vm);
 
    if ( vm.count("help") )
    {
      __ROOTCODE(
        std::cout << "Usage: " << argv[0] << " [options] file.xml\n"
                  << "  file.xml is an optional filename for XML input.\n"
                  << "  Default input is input.xml.\n\n"
                  << all << "\n"; 
      )
      return 1;
    }
    if ( vm.count("version") )
    {
      __ROOTCODE( 
        std::cout << "\n" << __PROGNAME__
                  << " from the " << PACKAGE_STRING << " package\n"
                  << "Subversion Revision: " << SVN::Revision << "\n\n"; 
      )
      return 1;
    }
    if ( vm.count("input") )
    {
      filename = vm["input"].as< std::string >();
      __ROOTCODE(std::cout << "Input: " << filename << ".\n";)
    }
    if ( vm.count("reruns") and vm["reruns"].as<unsigned>() > 1 )
    {
      reruns = vm["reruns"].as<unsigned>();
      __ROOTCODE( 
        std::cout << "Will rerun algorithm " << reruns 
                  << " independent times.\n";
      )
    }
      
    __MPICODE( boost::mpi::broadcast( ::mpi::main, filename, 0 ); )


    for( types::t_int i = 0; i < reruns; ++i )
    {
      GA::Darwin< t_Evaluator > ga;
      __DOASSERT( not ga.Load(filename),
                  "Could not load input from file " << filename << "\n" )
      
      std::cout << "Rerun " << i+1 << " of " << reruns << ".\n";
      ga.run();
      // for reruns.
      Print::xmg.dont_truncate();
      Print::out.dont_truncate();
    }
  }
  catch ( boost::program_options::invalid_command_line_syntax &_b)
  {
    std::cerr << "Caught error while running " << __PROGNAME__ << "\n"
              << "Something wrong with the command-line input.\n"
              << _b.what() << "\n";
    return 0;
  }
  catch ( boost::program_options::invalid_option_value &_i )
  {
    std::cerr << "Caught error while running " << __PROGNAME__ << "\n"
              << "Argument of option in command-line is invalid.\n"
              << _i.what() << "\n";
    return 0;
  }
  catch ( boost::program_options::unknown_option &_u)
  {
    std::cerr << "Caught error while running " << __PROGNAME__ << "\n"
              << "Unknown option in command-line.\n"
              << _u.what() << "\n";
    return 0;
  }
  catch (  boost::program_options::too_many_positional_options_error &_e )
  {
    std::cerr << "Caught error while running " << __PROGNAME__ << "\n"
              << "Too many arguments in command-line.\n"
              << _e.what() << "\n";
    return 0;
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
      std::string message;
      boost::mpi::gather( ::mpi::main, sstr.str(), message, 0 );
    )
    __NOTMPIROOT( return 0; )
    std::cerr << sstr.str() << std::endl;
    return 0;
  }
  return 1;
}

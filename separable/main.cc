//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <fstream>
#include <sstream>
#include <string>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/program_options.hpp>

#include <tinyxml/tinyxml.h> 

#include <opt/gsl_lsq.h>
#include "lsq.h"
#include "cefitting.h"

#include <revision.h>
#define __PROGNAME__ "Fixed-Lattice Sum of Separable function.\n" 

int main(int argc, char *argv[]) 
{
  try
  {
    namespace po = boost::program_options;

    po::options_description generic("Generic Options");
    generic.add_options()
           ("help,h", "produces this help message.")
           ("version,v", "prints version string.");
           ("verbose,p", po::value<bool>()->default_value(false),
                         "Verbose output.\n"  );
           ("reruns,r", po::value<unsigned>()->default_value(1),
                        "number of times to run the algorithm.\n" 
                        "Is equivalent to manually re-launching the program.\n");
    po::options_description specific("GA Options");
    specific.add_options()
        ("cross,c", po::value<types::t_real>()->default_value(1e-4),
                    "Performs leave-one-out"
                    " cross-validation, rather than simple fit.\n"  )
        ("size,s", po::value<types::t_unsigned>()->default_value(3),
                   "Size of the cubic basis." )
        ("rank,r", po::value<types::t_unsigned>()->default_value(3),
                   "Rank of the sum of separable functions." )
        ("tolerance", po::value<types::t_real>()->default_value(1e-4),
                      "Tolerance of the alternating linear-least square fit.\n"  )
        ("maxiter", po::value<types::t_unsigned>()->default_value(40),
                    "Maximum number of iterations for"
                    " Alternating linear-least square fit.\n"  )
        ("1dtolerance", po::value<types::t_real>()->default_value(1e-4),
                        "Tolerance of the 1d linear-least square fit.\n" )
        ("update", po::value<bool>()->default_value(false),
                   "Whether to update during 1d least-square fits.\n" )
        ("svd", po::value<bool>()->default_value(false),
                "Whether the 1d least-square fit should use"
                " single-value decomposition.\n"  );
 
    po::options_description all;
    all.add(generic).add(specific);
 
    po::positional_options_description p;
    p.add("dir", 1);
    p.add("input", 2);
 
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(all).positional(p).run(), vm);
    po::notify(vm);
 
    std::cout << "\n" << __PROGNAME__
              << " from the " << PACKAGE_STRING << " package\n"
              << "Subversion Revision: " << SVN::Revision << "\n\n"; 
    if ( vm.count("version") ) return 1;
    if ( vm.count("help") )
    {
      std::cout << "Usage: " << argv[0] << " [options] directory filename\n"
                  << "  directory is an optional path to the training set input.\n"
                  << "  Default is current directory.\n"
                  << "  filename is an option filename for the file from\n" 
                  << "  which to load the lattice.\n"
                  << "  (starting from current directory)\n"
                  << "  Default is input.xml.\n\n"
                  << all << "\n"; 
      return 1;
    }

    std::string dir(".");
    if( vm.count("dir") ) dir = vm["dir"].as< std::string >();
    std::string filename(".");
    if( vm.count("input") ) filename = vm["input"].as< std::string >();

    bool verbose = vm.count("verbose");
    types::t_unsigned reruns( vm["reruns"].as< types::t_unsigned >() );
    __ASSERT( reruns == 0, "0 number of runs performed... As required on input.\n" )
    bool cross = vm.count("cross");
    types::t_unsigned rank( vm["rank"].as< types::t_unsigned >() );
    types::t_unsigned size( vm["size"].as< types::t_unsigned >() );
    types::t_real tolerance( vm["tolerance"].as< types::t_real >() );
    types::t_unsigned maxiter( vm["maxiter"].as< types::t_unsigned >() );
    types::t_real dtolerance( vm["1dtolerance"].as< types::t_real >() );
    bool doupdate = vm.count("doupdate");
    bool svd = vm.count("svd");
    if( dtolerance > tolerance )
    {
      std::cout << "1d tolerance cannot be smaller then alternating tolerance.\n" 
                << tolerance << " < " << dtolerance << "\n";
      return 0;
    }

    std::cout << "Performing " << (cross ? "Cross-Validation" : "Fitting" ) << ".\n"
              << "Size of the cubic basis: " << size << "\n"
              << "Rank of the sum of separable functions: " << rank << "\n"
              << "Data directory: " << dir << "\n";
    if( reruns <= 1 )  std::cout << "single";
    else               std::cout << reruns;
    std::cout << " Run" << (reruns <= 1 ? ".": "s." )  << "\n"
              << (verbose ? "Verbose": "Quiet") << " output.\n"
              << "Alternating linear-least square tolerance: " 
                 << tolerance << "\n"
              << "Maximum number of iterations for alternating least-square fit: "
                 << maxiter << "\n"
              << "1d linear-least square tolerance: " 
                 << dtolerance << "\n"
              << "Will" << ( doupdate ? " ": " not " )
                 << "update between dimensions.\n"
              << ( svd ? "Single": "No Single" ) << " Value Decomposition."
              << std::endl;

    // Loads lattice
    { 
      TiXmlDocument doc( filename );
      __DOASSERT( not doc.LoadFile(), 
                    "error while opening lattice input file "
                 << filename << "\n" << doc.ErrorDesc()  )
      TiXmlHandle handle( &doc );
      TiXmlElement *child = handle.FirstChild( "Job" ).FirstChild( "Lattice" ).Element();
      __DOASSERT( not child, "Could not find Lattice in input." )
      Crystal::Structure::lattice = new Crystal :: Lattice;
      __DOASSERT( not Crystal::Structure::lattice->Load(*child),
                  "Error while reading Lattice from input.")
    }

    // Initializes fitting.
    Fitting::Allsq< Fitting::Gsl > allsq;
    allsq.itermax = maxiter;
    allsq.tolerance = tolerance;
    allsq.llsq.tolerance = dtolerance;
    allsq.llsq.dosvd = svd;
    allsq.llsq.doweights = true;

    // Initializes the symmetry-less separable function.
    CE::Separables separables( rank, size, "cube" );
    
    // Initializes cum-symmetry separable function.
    CE::SymSeparables symsep( separables );

    // Initializes collapse functor.
    Separable::Collapse< CE::Separables > collapse( separables );
    collapse.do_update = doupdate;

    // Initializes Interface to allsq.
    Fitting::SepCeInterface interface;
    interface.read( symsep, dir );
    
    // loops over reruns
    typedef Fitting::SepCeInterface::t_PairErrors t_PairErrors;
    if( reruns > 1 ) std::cout << "\n\n*********  Run 1.  *********\n";
    for( types::t_unsigned n(0); n < reruns; ++n )
    {
      if( n > 1 ) std::cout << "\n\n*********  Run " << n+1 << ".  *********\n";

      // fitting.
      if( not cross )
      {
        std::cout << "\nFitting using whole training set:" << std::endl;
        interface.fit( allsq, collapse );
        t_PairErrors result; 
        result = interface.check_training( separables, verbose );
        std::cout << " average error: " << result.first
                  << " maximum error: " << result.second << std::endl;
      }
      else
      {
        std::cout << "\nLeave-one-out prediction:" << std::endl;
        std::pair< t_PairErrors, t_PairErrors> result;
        result = Fitting::leave_one_out( interface, allsq, collapse, verbose );
        std::cout << " Training errors:\n"
                  << "    average error: " << result.first.first
                  << "    maximum error: " << result.first.second << "\n"
                  << " Prediction errors:\n"
                  << "    average error: " << result.second.first
                  << "    maximum error: " << result.second.second << std::endl;
      }
    }
    std::cout << "\n\n\nEnd of " << __PROGNAME__ << ".\n" << std::endl;

    if( Crystal::Structure::lattice ) delete Crystal::Structure::lattice;
    Crystal::Structure::lattice = NULL;
  }
  catch ( boost::program_options::invalid_command_line_syntax &_b)
  {
    if( Crystal::Structure::lattice ) delete Crystal::Structure::lattice;
    Crystal::Structure::lattice = NULL;
    std::cerr << "Caught error while running " << __PROGNAME__ << "\n"
              << "Something wrong with the command-line input.\n"
              << _b.what() << "\n";
    return 0;
  }
  catch ( boost::program_options::invalid_option_value &_i )
  {
    if( Crystal::Structure::lattice ) delete Crystal::Structure::lattice;
    Crystal::Structure::lattice = NULL;
    std::cerr << "Caught error while running " << __PROGNAME__ << "\n"
              << "Argument of option in command-line is invalid.\n"
              << _i.what() << "\n";
    return 0;
  }
  catch ( boost::program_options::unknown_option &_u)
  {
    if( Crystal::Structure::lattice ) delete Crystal::Structure::lattice;
    Crystal::Structure::lattice = NULL;
    std::cerr << "Caught error while running " << __PROGNAME__ << "\n"
              << "Unknown option in command-line.\n"
              << _u.what() << "\n";
    return 0;
  }
  catch (  boost::program_options::too_many_positional_options_error &_e )
  {
    if( Crystal::Structure::lattice ) delete Crystal::Structure::lattice;
    Crystal::Structure::lattice = NULL;
    std::cerr << "Caught error while running " << __PROGNAME__ << "\n"
              << "Too many arguments in command-line.\n"
              << _e.what() << "\n";
    return 0;
  }
  catch ( std::exception &e )
  {
    if( Crystal::Structure::lattice ) delete Crystal::Structure::lattice;
    Crystal::Structure::lattice = NULL;
    std::cerr << "Caught error while running " << __PROGNAME__ 
              << "\n" << e.what() << "\n";
    return 0;
  }
  return 1;
}


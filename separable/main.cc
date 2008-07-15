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
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/shared_ptr.hpp>

#include <tinyxml/tinyxml.h> 

#include <opt/allsq.h>
#include <opt/cgs.h>
#include <opt/types.h>
#include <opt/fuzzy.h>
#include <opt/debug.h>

#include "cefitting.h"
#include "cecollapse.h"
#include "bestof.h"
#include "leave_many_out.h"

#include <revision.h>
#define __PROGNAME__ "Fixed-Lattice Sum of Separable functions" 

const types::t_unsigned print_reruns   = 1;
const types::t_unsigned print_checks   = 2;
const types::t_unsigned print_allsq    = 3;
const types::t_unsigned print_function = 4;
const types::t_unsigned print_data     = 5;
const types::t_unsigned print_llsq     = 6;

int main(int argc, char *argv[]) 
{
  try
  {
    namespace po = boost::program_options;

    Fitting::LeaveManyOut leavemanyout;

    po::options_description generic("Generic Options");
    generic.add_options()
           ("help,h", "produces this help message.")
           ("version,v", "prints version string.")
           ("verbose,p", po::value<types::t_unsigned>()->default_value(0),
                         "Level of verbosity.\n"  )
           ("seed", po::value<types::t_unsigned>()->default_value(0),
                    "Seed of the random number generator.\n"  )
           ("reruns", po::value<types::t_unsigned>()->default_value(1),
                      "number of times to run the algorithm.\n" 
                      "Is equivalent to manually re-launching the program.\n");
    po::options_description specific("Separables Options");
    specific.add_options()
        ("cross,c", "Performs leave-one-out"
                    " cross-validation, rather than simple fit.\n"  )
        ("size,s", po::value<types::t_unsigned>()->default_value(3),
                   "Size of the cubic basis." )
        ("rank,r", po::value<types::t_unsigned>()->default_value(3),
                   "Rank of the sum of separable functions." )
        ("tolerance", po::value<types::t_real>()->default_value(1e-4),
                      "Tolerance of the alternating linear-least square fit.\n"  )
        ("maxiter,m", po::value<types::t_unsigned>()->default_value(40),
                      "Maximum number of iterations for"
                      " Alternating linear-least square fit.\n"  )
        ("1dtolerance", po::value<types::t_real>()->default_value(1e-4),
                        "Tolerance of the 1d linear-least square fit.\n" )
        ("noupdate", "Whether to update during 1d least-square fits.\n" )
        ("conv", "Use conventional cell rather than unit-cell.\n"
                 "Should work for fcc and bcc if lattice is inputed right.\n" )
        ("random", po::value<types::t_real>()->default_value(5e-1),
                   "Coefficients will be chosen randomly in the range [random,-random].\n" )
        ("nbguesses", po::value<types::t_unsigned>()->default_value(1),
                      "Number of initial guesses to try prior to (any) fitting.\n" );
    leavemanyout.add_cmdl( specific );
    po::options_description hidden("hidden");
    hidden.add_options()
        ("offset", po::value<types::t_real>()->default_value(0e0), 
                   "Adds an offset to the energies.\n" )
        ("prerun", "Wether to perform real-runs, or small pre-runs followed"
                   " by a a longer, converged run.\n" )
        ("datadir", po::value<std::string>()->default_value("./"))
        ("latticeinput", po::value<std::string>()->default_value("input.xml"));
 
    po::options_description all;
    all.add(generic).add(specific);
    po::options_description allnhidden;
    allnhidden.add(all).add(hidden);
 
    po::positional_options_description p;
    p.add("datadir", 1);
    p.add("latticeinput", 2);
 
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(allnhidden).positional(p).run(), vm);
    po::notify(vm);
 
    std::cout << "\n" << __PROGNAME__
              << " from the " << PACKAGE_STRING << " package.\n"
              << "Subversion Revision: " << SVN::Revision << "\n\n"; 
    if ( vm.count("version") ) return 1;
    if ( vm.count("help") )
    {
      std::cout << "Usage: " << argv[0] << " [options] DATADIR LATTICEINPUT\n"
                   "  _ DATATDIR (=./) is an optional path to the"
                   " training set input.\n"
                   "  _ LATTICEINPUT (=input.xml) is an optional filename"
                   " for the file\n"
                   "                 from which to load the lattice. "
                   "LATTICEINPUT should be\n"
                   "                 full path, a relative path "
                   "starting from the current\n"
                   "                 directory, or a relative path "
                   "starting from the DATADIR\n"
                   "                 directory (checked in that order.)\n\n" 
                << all << "\n"; 
      return 1;
    }

    std::string dir(".");
    if( vm.count("datadir") ) dir = vm["datadir"].as< std::string >();
    std::string filename("input.xml");
    if( vm.count("latticeinput") ) filename = vm["latticeinput"].as< std::string >();

    types::t_unsigned verbose = vm["verbose"].as<types::t_unsigned>();
    types::t_unsigned seed = vm["seed"].as<types::t_unsigned>();
    seed = opt::random::seed( seed );
    types::t_unsigned reruns(1);
    if( vm.count("reruns") ) reruns = vm["reruns"].as< types::t_unsigned >();
    __DOASSERT( reruns == 0, "0 number of runs performed... As required on input.\n" )
    bool cross = vm.count("cross");
    types::t_unsigned rank( vm["rank"].as< types::t_unsigned >() );
    __DOASSERT( rank == 0, "Separable function of rank 0 is obnoxious.\n" )
    types::t_unsigned size( vm["size"].as< types::t_unsigned >() );
    __DOASSERT( size == 0, "Separable function of dimension 0 is obnoxious.\n" )
    types::t_real tolerance( vm["tolerance"].as< types::t_real >() );
    types::t_unsigned maxiter( vm["maxiter"].as< types::t_unsigned >() );
    types::t_real dtolerance( vm["1dtolerance"].as< types::t_real >() );
    bool doupdate = not vm.count("noupdate");
    bool convcell = vm.count("conv");
    types::t_real offset ( vm["offset"].as< types::t_real >() );
    if( Fuzzy::eq( offset, types::t_real(0) ) ) offset = types::t_real(0);
    bool prerun ( vm.count("prerun") != 0 );
    types::t_real howrandom( vm["random"].as<types::t_real>() );
    types::t_unsigned nbguesses( vm["nbguesses"].as<types::t_unsigned>() );
    __ASSERT( nbguesses == 0, "Invalid input nbguesses = 0.\n" )

    // Loads lattice
    boost::shared_ptr< Crystal::Lattice > lattice( new Crystal :: Lattice );
    { 
      TiXmlDocument doc;
      if( boost::filesystem::exists( filename ) )
      {
        __DOASSERT( not doc.LoadFile( filename ), 
                     "Found " << filename << " but could not parse.\n"
                  << "Possible incorrect XML syntax.\n" 
                  << doc.ErrorDesc()  )
      }
      else 
      {
        boost::filesystem::path path( dir );
        boost::filesystem::path fullpath = path / filename;
        __DOASSERT( not boost::filesystem::exists( fullpath ),
                     "Could not find "<< filename 
                  << " in current directory, nor in " <<  path )
        __DOASSERT( not doc.LoadFile( fullpath.string() ),
                     "Could not parse " << fullpath 
                  << ".\nPossible incorrect XML syntax.\n"
                  << doc.ErrorDesc()  )
      }
      TiXmlHandle handle( &doc );
      TiXmlElement *child = handle.FirstChild( "Job" )
                                  .FirstChild( "Lattice" ).Element();
      __DOASSERT( not child, "Could not find Lattice in input." )
      Crystal::Structure::lattice = &( *lattice );
      __DOASSERT( not lattice->Load(*child),
                  "Error while reading Lattice from input.")
#     if defined (_TETRAGONAL_CE_)
        // Only Constituent-Strain expects and space group determination
        // expect explicitely tetragonal lattice. 
        // Other expect a "cubic" lattice wich is implicitely tetragonal...
        // Historical bullshit from input structure files @ nrel.
        for( types::t_int i=0; i < 3; ++i ) 
          if( Fuzzy::eq( lattice->cell.x[i][2], 0.5e0 ) )
            lattice->cell.x[i][2] = 0.6e0;
#     endif
      lattice->find_space_group();
#     if defined (_TETRAGONAL_CE_)
        // Only Constituent-Strain expects and space group determination
        // expect explicitely tetragonal lattice. 
        // Other expect a "cubic" lattice wich is implicitely tetragonal...
        // Historical bullshit from input structure files @ nrel.
        for( types::t_int i=0; i < 3; ++i ) 
          if( Fuzzy::eq( lattice->cell.x[i][2], 0.6e0 ) )
            lattice->cell.x[i][2] = 0.5e0;
#     endif
    }

    // Initializes fitting.
    typedef Fitting::Allsq<Fitting::Cgs> t_Fitting;
    t_Fitting allsq;

    Separables::BestOf< t_Fitting > bestof;
    bestof.n = reruns;
    bestof.verbose = verbose >= print_reruns;
    bestof.prerun = prerun;

    bestof.t_Fitting::itermax = maxiter;
    bestof.t_Fitting::tolerance = tolerance;
    bestof.t_Fitting::verbose = verbose >= print_allsq;
    bestof.t_Fitting::do_update = doupdate;
    bestof.t_Fitting::linear_solver.tolerance = dtolerance;
    bestof.t_Fitting::linear_solver.verbose = verbose >= print_llsq;

    // Initializes the symmetry-less separable function.
    typedef CE::Separables t_Function;
    t_Function separables( rank, size, convcell ? "conv": "cube" );
    
    // Initializes cum-symmetry separable function.
    CE::SymSeparables symsep( separables );

    // Initializes collapse functor.
    Separable::EquivCollapse< t_Function > collapse( separables );

    // Initializes Interface to allsq.
    Fitting::SepCeInterface interface;
    interface.howrandom = howrandom;
    interface.nb_initial_guesses = nbguesses;
    interface.verbose = verbose >= 4;
    interface.set_offset( offset );
    interface.read( symsep, dir, "LDAs.dat", verbose >= print_data );
#   if defined (_TETRAGONAL_CE_)
      // From here on, lattice should be explicitely tetragonal.
      for( types::t_int i=0; i < 3; ++i ) 
        if( Fuzzy::eq( lattice->cell.x[i][2], 0.5e0 ) )
          lattice->cell.x[i][2] = 0.6e0;
#   endif
    types::t_real mean = interface.mean();
    types::t_real variance = interface.variance();

    // extract leave-many-out commandline
    leavemanyout.extract_cmdl( vm );
    leavemanyout.verbose = verbose;

    

    typedef Fitting::SepCeInterface::t_PairErrors t_PairErrors;

    std::cout << "Performing " << (cross ? "Cross-Validation" : "Fitting" ) << ".\n"
              << "Using " << ( convcell ? "conventional ": "unit-" )
                 << "cell for basis determination.\n"
              << "Size of a separable function " << separables.size() << "\n"
              << "Rank of the sum of separable functions: " << rank << "\n"
              << "d.o.f.: " << separables.size() * rank << "\n"
              << "Data directory: " << dir << "\n";
    if( reruns <= 1 )  std::cout << "single";
    else               std::cout << reruns;
    std::cout << " Run" << (reruns <= 1 ? ".": "s." )  << "\n";
    if( not verbose ) std::cout << "Quiet output.\n";
    else std::cout << "Level of verbosity: " << verbose << "\n";
    std::cout << "Alternating linear-least square tolerance: " 
                 << tolerance << "\n"
              << "Maximum number of iterations for alternating least-square fit: "
                 << maxiter << "\n"
              << "Alternating linear-least square tolerance: " 
                 << tolerance << "\n"
              << "1d linear-least square tolerance: " 
                 << dtolerance << "\n"
              << "Will" << ( doupdate ? " ": " not " )
                 << "update between dimensions.\n"
              << "Data mean: " << mean << "\n"
              << "Data Variance: " << variance << "\n"
              << "Range of initial guesses:[ " <<  howrandom << ", " << howrandom << "].\n"
              << "Number of initial guesses: " <<  nbguesses << ".\n";
    if( prerun )
     std::cout << "Performing prerun.\n";
    std::cout << "Random Seed: " << seed << "\n";

    if( Fuzzy :: neq( offset, 0e0 ) ) std::cout << "Offset: " << offset << "\n";

    // fitting.
    if( leavemanyout.do_perform )
    {
      std::cout << "\nStarting leave-many out predictive fit." << std::endl;
      Fitting::LeaveManyOut::t_Return result;
      result = leavemanyout( interface, bestof, collapse );
      std::cout << " Training errors:\n"
                << "    average error: " << result.first.first
                <<  " ( " << result.first.first / std::abs(mean) *1e2 << "% )" 
                << "    maximum error: " << result.first.second << "\n"
                <<  " ( " << result.first.second / std::abs(mean) *1e2 << "% )\b" 
                << " Prediction errors:\n"
                << "    average error: " << result.second.first
                <<  " ( " << result.second.first / std::abs(mean) *1e2 << "% )" 
                << "    maximum error: " << result.second.second 
                <<  " ( " << result.second.second / std::abs(mean) *1e2 << "% )\n"; 
    }
    else if( not cross )
    {
      std::cout << "\nFitting using whole training set:" << std::endl;
      types::t_real residual = interface.fit( bestof, collapse );
      t_PairErrors result; 
      result = interface.check_training( separables, verbose >= print_checks );
      std::cout << " Residual: " << residual 
                <<  " ( " << residual / variance *1e2 << "% )" 
                << " average error: " << result.first
                <<  " ( " << result.first / std::abs(mean) *1e2 << "% )" 
                << " maximum error: " << result.second 
                <<  " ( " << result.second / std::abs(mean) *1e2 << "% )" 
                << std::endl;
    }
    else
    {
      std::cout << "\nLeave-one-out prediction:" << std::endl;
      std::pair< t_PairErrors, t_PairErrors> result;
      result = Fitting::leave_one_out( interface, bestof, collapse, verbose >= print_checks );
      std::cout << " Training errors:\n"
                << "    average error: " << result.first.first
                <<  " ( " << result.first.first / std::abs(mean) *1e2 << "% )" 
                << "    maximum error: " << result.first.second << "\n"
                <<  " ( " << result.first.second / std::abs(mean) *1e2 << "% )\n" 
                << " Prediction errors:\n"
                << "    average error: " << result.second.first
                <<  " ( " << result.second.first / std::abs(mean) *1e2 << "% )" 
                << "    maximum error: " << result.second.second
                <<  " ( " << result.second.second / std::abs(mean) *1e2 << "% )\n"; 
    }
    if( verbose >= print_function ) std::cout << separables << "\n";

    std::cout << "\n\n\nEnd of " << __PROGNAME__ << ".\n" << std::endl;

  }
  catch ( boost::program_options::invalid_command_line_syntax &_b)
  {
    std::cerr << "Caught error while running " << __PROGNAME__ << "\n"
              << "Something wrong with the command-line input.\n"
              << _b.what() << std::endl;
    return 0;
  }
  catch ( boost::program_options::invalid_option_value &_i )
  {
    std::cerr << "Caught error while running " << __PROGNAME__ << "\n"
              << "Argument of option in command-line is invalid.\n"
              << _i.what() << std::endl;
    return 0;
  }
  catch ( boost::program_options::unknown_option &_u)
  {
    std::cerr << "Caught error while running " << __PROGNAME__ << "\n"
              << "Unknown option in command-line.\n"
              << _u.what() << std::endl;
    return 0;
  }
  catch (  boost::program_options::too_many_positional_options_error &_e )
  {
    std::cerr << "Caught error while running " << __PROGNAME__ << "\n"
              << "Too many arguments in command-line.\n"
              << _e.what() << std::endl;
    return 0;
  }
  catch ( std::exception &e )
  {
    std::cerr << "Caught error while running " << __PROGNAME__ 
              << "\n" << e.what() << std::endl;
    return 0;
  }
  return 1;
}


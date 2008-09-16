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
#include <boost/bind.hpp>

#include <opt/cgs.h>
#include <opt/types.h>
#include <opt/fuzzy.h>
#include <opt/debug.h>
#include <opt/errors.h>
#include <opt/random.h>
#include <crystal/lattice.h>
#include <crystal/structure.h>

#include "functional.h"
#include "sepmappings.h"
#include "colmappings.h"
#include "collapse.h"
#include "prepare.h"
#include "methods.h"

#include <revision.h>
#define __PROGNAME__ "Fixed-Lattice Sum of Separable functions" 

const types::t_unsigned print_reruns   = 1;
const types::t_unsigned print_checks   = 2;
const types::t_unsigned print_function = 3;
const types::t_unsigned print_allsq    = 4;
const types::t_unsigned print_data     = 5;
const types::t_unsigned print_llsq     = 6;

int main(int argc, char *argv[]) 
{
  try
  {
    namespace po = boost::program_options;
    namespace bl = boost::lambda;
    namespace fs = boost::filesystem;

    po::options_description generic("Generic Options");
    generic.add_options()
           ("help,h", "produces this help message.")
           ("version,v", "prints version string.")
           ("verbose,p", po::value<types::t_unsigned>()->default_value(0),
                         "Level of verbosity.\n"  )
           ("seed", po::value<types::t_unsigned>()->default_value(0),
                    "Seed of the random number generator.\n"  );
    po::options_description specific("Separables Options");
    specific.add_options()
        ("loo,l", "Performs leave-one-out cross-validation.\n"  )
        ("nofit,l", "Does not perform fit.\n"  )
        ("rank,r", po::value<types::t_unsigned>()->default_value(3),
                   "Rank of the sum of separable functions." )
        ("basis,r", po::value<std::string>()->default_value("1x1x4"),
                   "Description of the ranks/size of the figure used\n." )
        ("tolerance", po::value<types::t_real>()->default_value(1e-4),
                      "Tolerance of the alternating linear-least square fit.\n"  )
        ("maxiter,m", po::value<types::t_unsigned>()->default_value(40),
                      "Maximum number of iterations for"
                      " Alternating linear-least square fit.\n"  )
        ("1dtolerance", po::value<types::t_real>()->default_value(1e-4),
                        "Tolerance of the 1d linear-least square fit.\n" ) 
        ("random", po::value<types::t_real>()->default_value(5e-1),
                   "Coefficients' randomness.\n" )
        ("lambda,l", po::value<types::t_real>()->default_value(0),
                     "Regularization factor.\n" )
        ("bestof,b", po::value<types::t_unsigned>()->default_value(1),
                     "Performs best-of fit.\n" )
        ("which,w", po::value<types::t_unsigned>()->default_value(0),
                     "Performs best-of for 0 (variance), 1 (mean), or 2(max).\n" );
    po::options_description hidden("hidden");
    hidden.add_options()
        ("datadir", po::value<std::string>()->default_value("./"))
        ("latinput", po::value<std::string>()->default_value("input.xml"));
 
    po::options_description all;
    all.add(generic).add(specific);
    po::options_description allnhidden;
    allnhidden.add(all).add(hidden);
 
    po::positional_options_description p;
    p.add("datadir", 1);
    p.add("latinput", 2);
 
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
                   "                 a full path or a relative path "
                   "starting from the current\n"
                   "                 directory, or a relative path "
                   "starting from the DATADIR\n"
                   "                 directory (checked in that order.)\n\n" 
                << all << "\n"; 
      return 1;
    }

    fs::path dir( vm["datadir"].as< std::string >() );
    __DOASSERT( not fs::exists( dir ), dir << " does not exist.\n" );
    __DOASSERT( not fs::exists( dir / "LDAs.dat" ), 
                ( dir / "LDAs.dat" ) << " does not exist.\n" );
    __DOASSERT
    (
      not (    fs::is_regular( dir / "LDAs.dat" ) 
            or fs::is_symlink( dir / "LDAs.dat" ) ),
      ( dir / "LDAs.dat" ) << " is neither a regular file nor a symlink.\n"
    ) 
    fs::path latinput( vm["latinput"].as< std::string >() );
    if(    ( not fs::exists( latinput ) )
        or ( not ( fs::is_symlink(latinput) or fs::is_regular(latinput) ) ) )
     latinput = dir / latinput;
    __DOASSERT( not fs::exists( latinput ),
                   "Lattice input file " << latinput
                << " could not be found in ./ nor in " << dir << ".\n" )
    __DOASSERT( not ( fs::is_regular( latinput ) or fs::is_symlink( latinput ) ),
                latinput << " is a not a valid file.\n" );

    const bool dofit( vm.count("nofit") == 0 );
    const bool doloo( vm.count("loo") != 0 );
    __ASSERT( not ( dofit or doloo ), "Nothing to do..." )
    const types::t_unsigned verbosity = vm["verbose"].as<types::t_unsigned>();
    types::t_unsigned seed = vm["seed"].as<types::t_unsigned>();
    seed = opt::random::seed( seed );
    const types::t_unsigned rank( vm["rank"].as< types::t_unsigned >() );
    const std::string bdesc( vm["basis"].as<std::string>() );
    __DOASSERT( rank == 0, "Separable function of rank 0 is obnoxious.\n" )
    const types::t_real tolerance( vm["tolerance"].as< types::t_real >() );
    const types::t_unsigned maxiter( vm["maxiter"].as< types::t_unsigned >() );
    const types::t_real dtolerance( vm["1dtolerance"].as< types::t_real >() );
    const types::t_real howrandom( vm["random"].as<types::t_real>() );
    const types::t_real lambda( vm["lambda"].as<types::t_real>() );
    const types::t_unsigned bestof( vm["bestof"].as<types::t_unsigned>() );
    __DOASSERT( bestof == 0, "0 jobs to be performed..." )
    const types::t_unsigned which( vm["which"].as<types::t_unsigned>() );
    __DOASSERT( which >= 3, "Don't know which error to perform bestof for.\n" )

    // Loads lattice
    boost::shared_ptr< Crystal::Lattice >
    lattice( Crystal::read_lattice( latinput, dir ) );
    Crystal::Structure::lattice = lattice.get();
#   if defined (_TETRAGONAL_CE_)
      // Only Constituent-Strain expects and space group determination
      // expect explicitely tetragonal lattice. 
      // Other expect a "cubic" lattice wich is implicitely tetragonal...
      // Historical bullshit from input structure files @ nrel.
      for( types::t_int i=0; i < 3; ++i ) 
        if( Fuzzy::eq( lattice->cell.x[i][2], 0.5e0 ) )
          lattice->cell.x[i][2] = 0.6e0;
#   endif
    lattice->find_space_group();
#   if defined (_TETRAGONAL_CE_)
      // Only Constituent-Strain expects and space group determination
      // expect explicitely tetragonal lattice. 
      // Other expect a "cubic" lattice wich is implicitely tetragonal...
      // Historical bullshit from input structure files @ nrel.
      for( types::t_int i=0; i < 3; ++i ) 
        if( Fuzzy::eq( lattice->cell.x[i][2], 0.6e0 ) )
          lattice->cell.x[i][2] = 0.5e0;
#   endif

    // Reads structures.
    std::vector< Crystal::Structure > structures;
    Crystal::read_ce_structures( dir / "LDAs.dat", structures );
    if( verbosity >= print_data )
      std::for_each
      (
        structures.begin(), structures.end(), 
        std::cout << bl::_1 << "\n"
      );


    // Initializes fitting.
    typedef Fitting::AlternatingLeastSquare<Fitting::Cgs> t_Fitting;
    t_Fitting allsq;
    allsq.itermax = maxiter;
    allsq.tolerance = tolerance;
    allsq.verbose = verbosity >= print_allsq;
    allsq.linear_solver.tolerance = dtolerance;
    allsq.linear_solver.verbose = verbosity >= print_llsq;

    // Initializes basis.
    CE::PosToConfs postoconfs( *Crystal::Structure::lattice );
    postoconfs.create_positions( bdesc );

    // Initializes the symmetry-less separable function.
    typedef Traits::CE::Separables< CE::Mapping::VectorDiff<2> > t_FunctionTraits;
    typedef CE::Separables< t_FunctionTraits > t_Function;
    t_Function separables;
    separables.set_rank_n_size( rank, postoconfs.positions.size() );
    separables.randomize( howrandom );
    separables.normalize();

    // Initializes collapse functor.
    typedef CE::Mapping::ExcludeOne< CE::Mapping::SymEquiv > t_Mapping;
    typedef CE::Policy::Regularization< t_Function > t_Regularization;
    typedef boost::numeric::ublas::matrix<size_t> t_Confs;
    typedef CE::Policy::HighMemUpdate< t_Function, t_Mapping, t_Confs > t_UpdatePolicy;
    typedef Traits::CE::Collapse< t_Function, t_Mapping, 
                                  t_Regularization, t_Confs,
                                  t_UpdatePolicy > t_CollapseTraits;
    typedef CE::Collapse< t_CollapseTraits > t_Collapse;
    t_Collapse collapse;
    collapse.init( structures, postoconfs );
    collapse.regularization().lambda = lambda;
    collapse.init( separables );


    opt::NErrorTuple nerror( opt::mean_n_var(structures) ); 

    std::cout << "Shape of separable function: " << bdesc << "\n"
              << "Rank of separable function " << rank << "\n"
              << "Size of separable function "
              << postoconfs.positions.size() << "\n"
              << "Data directory: " << dir << "\n";
    if( not verbosity ) std::cout << "Quiet output.\n";
    else std::cout << "Level of verbosity: " << verbosity << "\n";
    std::cout << "Alternating linear-least square tolerance: " 
                 << tolerance << "\n"
              << "Maximum number of iterations for alternating least-square fit: "
                 << maxiter << "\n"
              << "1d linear-least square tolerance: " 
                 << dtolerance << "\n"
              << "Data mean: " << nerror.nmean() << "\n"
              << "Data Variance: " << nerror.nvariance() << "\n"
              << "Random Seed: " << seed << "\n"
              << "Regulation factor: " << lambda << "\n"
              << "Performing best of  " << bestof;
    if( which == 0 ) std::cout << " for variance\n";
    if( which == 1 ) std::cout << " for mean\n";
    if( which == 2 ) std::cout << " for max\n";


    // Initializes best of fit.
    CE::Method::Fit< CE::Method::Policy::BestOf< CE::CollapseState > > fit( structures );
    fit.verbosity = verbosity -1;
    fit.policy.restarts = bestof;
    fit.policy.which = which;

    // fitting.
    if( doloo )
    {
      collapse.mapping().do_exclude = true;
      std::cout << "Starting Leave-One-Out Procedure.\n";
      opt::t_ErrorPair errors;
      errors = CE::Method::leave_one_out( collapse, fit, allsq, verbosity - 1 );
      std::cout << "Average Training Errors:\n " << ( nerror = errors.first ) << "\n";
      std::cout << "Final Prediction Errors:\n " << ( nerror = errors.second ) << "\n\n";
      std::cout << "Average Training Errors:\n " << errors.first << "\n";
      std::cout << "Final Prediction Errors:\n " << errors.second << "\n\n";
      if( verbosity >= print_function ) std::cout << separables << "\n";
    }
    if( dofit )
    {
      collapse.mapping().do_exclude = false;
      std::cout << "\nFitting using whole training set:" << std::endl;
      nerror = fit( collapse, allsq );
      std::cout << nerror << "\n"; 
      if( verbosity >= print_function ) std::cout << separables << "\n";
    }
    std::cout << "\n\n\nEnd of " << __PROGNAME__ << ".\n" << std::endl;

  }
  catch ( boost::program_options::invalid_command_line_syntax &_b)
  {
    std::cout << "Caught error while running " << __PROGNAME__ << "\n"
              << "Something wrong with the command-line input.\n"
              << _b.what() << std::endl;
    return 0;
  }
  catch ( boost::program_options::invalid_option_value &_i )
  {
    std::cout << "Caught error while running " << __PROGNAME__ << "\n"
              << "Argument of option in command-line is invalid.\n"
              << _i.what() << std::endl;
    return 0;
  }
  catch ( boost::program_options::unknown_option &_u)
  {
    std::cout << "Caught error while running " << __PROGNAME__ << "\n"
              << "Unknown option in command-line.\n"
              << _u.what() << std::endl;
    return 0;
  }
  catch (  boost::program_options::too_many_positional_options_error &_e )
  {
    std::cout << "Caught error while running " << __PROGNAME__ << "\n"
              << "Too many arguments in command-line.\n"
              << _e.what() << std::endl;
    return 0;
  }
  catch ( std::exception &e )
  {
    std::cout << "Caught error while running " << __PROGNAME__ 
              << "\n" << e.what() << std::endl;
    return 0;
  }
  return 1;
}


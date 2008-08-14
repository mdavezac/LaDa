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

#include <opt/types.h>
#include <opt/debug.h>
#include <opt/random.h>
#include <opt/gsl_mins.h>
#include <crystal/lattice.h>
#include <crystal/structure.h>

#include <revision.h>
#define __PROGNAME__ "Regulated Figure Search."

#include "regularization.h"
#include "cluster.h"

const types::t_int least = 0;
const types::t_int serrors = 1;
const types::t_int detailederrors = 2;
const types::t_int outermin = 3;
const types::t_int dclusters = 5;
const types::t_int startclusters = 6;
const types::t_int startstructures = 6;
const types::t_int allclusters = 7;
const types::t_int innermin = 8;

int main(int argc, char *argv[]) 
{
  namespace fs = boost::filesystem;
  namespace bl = boost::lambda;
  try
  {
    namespace po = boost::program_options;

    po::options_description generic("Generic Options");
    generic.add_options()
      ("help,h", "produces this help message.")
      ("version,v", "prints version string.")
      ("verbose,p", po::value<types::t_int>()->default_value(0),
                    "Level of verbosity."  )
      ("seed", po::value<types::t_unsigned>()->default_value(0),
               "Seed of the random number generator."  );
    po::options_description specific("Regularization Options");
    specific.add_options()
      ("latinput", po::value<std::string>()->default_value("input.xml"),
                   "Lattice input file." )
      ("tolerance", po::value<types::t_real>()->default_value(1e-4),
                    "Tolerance of the non-linear-least square fit."  )
      ("itermax,i", po::value<types::t_unsigned>()->default_value(0),
                    "Maximum number of iterations for the minimizer."  )
      ("maxpairs,m", po::value<types::t_unsigned>()->default_value(5),
                     "Max distance for pairs (in neighbors)."  )
      ("minimizer", po::value<std::string>()->default_value("simplex"),
                    "Type of Minimizer"  )
      ("iw", po::value<types::t_real>()->default_value(0),
             "Initial value of the weights."  )
      ("loo", "Single leave-one-out procedure." )
      ("fit", "Single leave-none-out procedure." )
      ("alpha", po::value<types::t_real>()->default_value(1),
                 "Lambda for pair regularization. With --loo/fit and --tcoef only. " )
      ("tcoef", po::value<types::t_real>()->default_value(0),
                "\"t\" coefficient. With --loo/fit only. " )
      ("volkerreg", " Pair regulation done with Volker Blum's"
                    " normalization (changes the units of the \"t\" coefficient). " );
       
    po::options_description hidden("hidden");
    hidden.add_options()
        ("datadir", po::value<std::string>()->default_value("./"))
        ("jtypes", po::value<std::string>()->default_value("jtypes"));
 
    po::options_description all;
    all.add(generic).add(specific);
    po::options_description allnhidden;
    allnhidden.add(all).add(hidden);
 
    po::positional_options_description p;
    p.add("datadir", 1);
    p.add("jtypes", 2);
 
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
      std::cout << "Usage: regular  [options] DATADIR JTYPES\n"
                   "  _ DATATDIR (=./) is an optional path to the training set input.\n"
                   "  _ JTYPES (=jtypes) is an optional filename for the file"
                   "from which to load\n"
                   "                     the lattice. (DATADIR must be present"
                   "when using JTYPES).\n"
                   "                     JTYPES should be a full path or a"
                   "relative path\n"
                   "                     starting from the current directory,"
                   " or a relative path\n"
                   "                     starting from the DATADIR directory"
                   "(checked in that\n"
                   "                     order.)\n"
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

    fs::path jtypes( vm["jtypes"].as< std::string >() );
    if(    ( not fs::exists( jtypes ) )
        or ( not ( fs::is_symlink(jtypes) or fs::is_regular(jtypes) ) ) )
      jtypes = dir / jtypes;
    __DOASSERT( not fs::exists( jtypes ),
                   "Jtypes input file " << jtypes
                << " could not be found in ./ nor in " << dir << ".\n" )
    __DOASSERT( not ( fs::is_regular( jtypes ) or fs::is_symlink( jtypes ) ),
                jtypes << " is a not a valid file.\n" );

    fs::path latinput( vm["latinput"].as< std::string >() );
    if(    ( not fs::exists( latinput ) )
        or ( not ( fs::is_symlink(latinput) or fs::is_regular(latinput) ) ) )
     latinput = dir / latinput;
    __DOASSERT( not fs::exists( latinput ),
                   "Lattice input file " << latinput
                << " could not be found in ./ nor in " << dir << ".\n" )
    __DOASSERT( not ( fs::is_regular( latinput ) or fs::is_symlink( latinput ) ),
                latinput << " is a not a valid file.\n" );

    types::t_unsigned maxpairs = vm["maxpairs"].as<types::t_unsigned>();
    types::t_int verbosity = vm["verbose"].as<types::t_int>();
    types::t_unsigned seed = vm["seed"].as<types::t_unsigned>();
    seed = opt::random::seed( seed );
    types::t_real tolerance( vm["tolerance"].as< types::t_real >() );
    types::t_unsigned itermax( vm["itermax"].as< types::t_unsigned >() );
    std::string minimizer_type( vm["minimizer"].as< std::string >() );
    types::t_real iweights( vm["iw"].as< types::t_real >() );
    bool doloo = vm.count("loo") != 0;
    bool dofit = vm.count("fit") != 0;
    bool volkerreg = vm.count("volkerreg") != 0;
    std::pair<types::t_real, types::t_real> pairreg( vm["alpha"].as<types::t_real>(),
                                                     vm["tcoef"].as<types::t_real>() );

    std::cout << " Input Parameters:\n"
              << "   Verbosity: " << verbosity << "\n"
              << "   Seed: " << seed << "\n"
              << "   Tolerance: " << tolerance << "\n"
              << "   Maximum number of iterations: " << itermax << "\n"
              << "   Lattice input file: " << latinput << "\n"
              << "   LDA energy input file: LDAs.dat\n"
              << "   Directory of structures input files: " << dir << "\n"
              << "   J0, J1, and many-body description file: " << jtypes << "\n"
              << "   Maximum distance of pairs (nth neighbor): " << maxpairs << "\n";
    if( not ( doloo or dofit ) )
      std::cout << "   Non-linear Minimizer: " << minimizer_type << "\n"
                << "   Initial value of the weights: " << iweights << "\n";
    else
    {
      if( doloo ) std::cout << "   Will perform leave-one-out.\n";
      if( dofit ) std::cout << "   Will perform leave-none-out.\n";
      std::cout << "   alpha: " << pairreg.first << "\n"
                << "   t: " << pairreg.second << "\n";
      if( volkerreg ) std::cout << "   Using Volker's normalization for t.\n";
    }
    std::cout << "\n";

    // consistency checks.
    __ASSERT( Fuzzy::le(pairreg.second, 0e0), "Coefficient \"t\" cannot negative.\n" )

    // Loads lattice
    boost::shared_ptr< Crystal::Lattice >
      lattice( Crystal::read_lattice( latinput, dir ) );
    Crystal::Structure::lattice = lattice.get();

    // Construct regularization.
    CE::Regulated reg;
    reg.cgs.tolerance = tolerance;
    reg.cgs.verbose = verbosity >= innermin;
    reg.cgs.itermax = 40;
    
    // reads jtypes
    CE::Regulated :: t_Clusters clusters;
    // add pair terms.
    CE::create_pairs( *lattice, maxpairs, clusters );
    std::cout << "Creating " << clusters.size() << " pair figures.\n";
    CE::read_clusters( *lattice, jtypes, reg.clusters ); 
    std::cout << "Read " << reg.clusters.size()
              << " figures from " << jtypes << ".\n";
    reg.clusters.insert(   reg.clusters.begin()
                         + std::min( (size_t)2, reg.clusters.size() ), 
                         clusters.begin(), clusters.end() );
    clusters.clear();
    if( verbosity >= startclusters )
    {
      CE::Regulated :: t_Clusters :: const_iterator i_class = reg.clusters.begin();
      CE::Regulated :: t_Clusters :: const_iterator i_class_end = reg.clusters.end();
      for(; i_class != i_class_end; ++i_class )
      {
        if( verbosity < allclusters )  std::cout << i_class->front()
                                                 << " D=" << i_class->size() 
                                                 << "\n";
        else std::for_each( 
                            i_class->begin(), i_class->end(), 
                            std::cout << bl::_1 << "\n"
                          );
      }
    }

    // now for the real job.
    if( doloo or dofit )
    {
      // Create object for fitting procedures.
      CE::Fit< CE::FittingPolicy::PairReg<> > fit;
      fit.alpha = pairreg.first;
      fit.tcoef = pairreg.second;
      fit.do_pairreg = ( not Fuzzy::is_zero( pairreg.second ) and maxpairs );
      fit.laksreg = not volkerreg;
      fit.verbose = verbosity >= detailederrors;
      reg.cgs.verbose = verbosity >= outermin;
      // initialization
      Crystal::read_ce_structures( dir / "LDAs.dat", fit.structures );
      if( verbosity >= startstructures )
        std::for_each
        (
          fit.structures.begin(), fit.structures.end(), 
          std::cout << bl::_1 << "\n"
        );
      fit.init( reg.clusters );
      opt::NErrorTuple nerror = fit.mean_n_var();

      CE::BaseFit::t_Vector x( reg.clusters.size() );
      if( doloo )
      {
        std::cout << "Starting Leave-One-Out Procedure.\n";
        std::pair< opt::ErrorTuple, opt::ErrorTuple > errors;
        errors = leave_one_out( fit, reg.cgs, x, verbosity >= serrors );
        std::cout << "Average Training Errors:\n " << ( nerror = errors.first ) << "\n";
        std::cout << "Final Prediction Errors:\n " << ( nerror = errors.second ) << "\n\n";
      }
      if( dofit )
      {
        std::cout << "Starting Fitting Procedure.\n";
        opt::ErrorTuple errors( fit( x, reg.cgs ) );
        std::cout << "Training Errors:\n " << ( nerror = errors ) << "\n\n";
      }
    }
    else
    {
      // initialization.
      Crystal::read_ce_structures( dir / "LDAs.dat", reg.structures );
      if( verbosity >= startstructures )
        std::for_each
        (
          reg.structures.begin(), reg.structures.end(), 
          std::cout << bl::_1 << "\n"
        );
      reg.init();
      opt::NErrorTuple nerror = reg.mean_n_var();
    
      if( minimizer_type.compare( "simplex" ) == 0 )
      {
        Minimizer::Simplex simplex;
        simplex.tolerance = tolerance;
        simplex.verbose = verbosity >= outermin;
        simplex.itermax = itermax;
        simplex.stepsize = 1;
        CE::drautz_diaz_ortiz( reg, simplex, verbosity, iweights);
      }
      else 
      {
        Minimizer::Gsl gsl;
        gsl.type =  Minimizer::Gsl::SteepestDescent;
        gsl.tolerance = tolerance;
        gsl.verbose = verbosity >= outermin;
        gsl.itermax = itermax;
        gsl.linestep = 0.01;
        gsl.linetolerance = tolerance * 1e1;
        if( minimizer_type.compare("bfgs2") == 0 )
          gsl.type = Minimizer::Gsl::BFGS2;
        if( minimizer_type.compare("bfgs") == 0 )
          gsl.type = Minimizer::Gsl::BFGS;
        if( minimizer_type.compare("sd") == 0 )
          gsl.type = Minimizer::Gsl::SteepestDescent;
        if( minimizer_type.compare("fr") == 0 )
          gsl.type = Minimizer::Gsl::FletcherReeves;
        if( minimizer_type.compare("pr") == 0 )
          gsl.type = Minimizer::Gsl::PolakRibiere;
      
        CE::drautz_diaz_ortiz( reg, gsl, verbosity, iweights);
      }
    }
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


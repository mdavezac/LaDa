#include "LaDaConfig.h"

#include <fstream>
#include <sstream>
#include <string>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/program_options.hpp>
#include <boost/mpl/vector.hpp>
#ifdef LADA_MPI
# include <boost/scoped_ptr.hpp>
# include <boost/mpi/environment.hpp>
#endif

# include <opt/mpi.h>
# include <opt/mpi.h>
#include <opt/types.h>
#include <opt/debug.h>
#include <math/random.h>
#include <opt/bpo_macros.h>
#include <minimizer/frprmn.h>
#ifdef LADA_WITH_GSL
# include <minimizer/gsl_mins.h>
#endif
#ifdef LADA_WITH_MINUIT2
# include <minimizer/minuit2.h>
#endif
#if defined(LADA_WITH_GSL) or defined(LADA_WITH_MINUIT2)
# include <minimizer/variant.h>
#endif
#include <crystal/lattice.h>
#include <crystal/structure.h>
#include <crystal/read_structure.h>

#define __PROGNAME__ "Regulated Figure Search."

#include "drautz_diaz_ortiz.h"
#include "regularization.h"
#include "cluster.h"
#include "create_pairs.h"

const LaDa::types::t_int least = 0;
const LaDa::types::t_int printinteractionenergies = 1;
const LaDa::types::t_int serrors = 1;
const LaDa::types::t_int detailederrors = 2;
const LaDa::types::t_int outermin = 3;
const LaDa::types::t_int dclusters = 5;
const LaDa::types::t_int startclusters = 6;
const LaDa::types::t_int startstructures = 6;
const LaDa::types::t_int allclusters = 7;
const LaDa::types::t_int innermin = 8;

int main(int argc, char *argv[]) 
{
  namespace fs = boost::filesystem;
  namespace bl = boost::lambda;
  namespace po = boost::program_options;
  LADA_TRY_BEGIN
  LADA_MPI_START


  LaDa::Fitting::LeaveManyOut leavemanyout;

  __BPO_START__
      ("verbose,p", po::value<LaDa::types::t_int>()->default_value(0),
                    "Level of verbosity."  )
      ("seed", po::value<LaDa::types::t_unsigned>()->default_value(0),
               "Seed of the random number generator."  );
  __BPO_SPECIFICS__("Regularization Options")
      ("latinput", po::value<std::string>()->default_value("input.xml"),
                   "Lattice input file." )
      ("tolerance", po::value<LaDa::types::t_real>()->default_value(1e-4),
                    "Tolerance of the non-linear-least square fit."  )
      ("linetol", po::value<LaDa::types::t_real>()->default_value(1e-4),
                    "Line tolerance for \"some\" minimizers." )
      ("linestep", po::value<LaDa::types::t_real>()->default_value(1e-4),
                    "Line step for \"some\" minimizers." )
      ("strategy", po::value<std::string>()->default_value("fast"),
                    "Strategy for Minuit2 minimizer." )
      ("itermax,i", po::value<LaDa::types::t_unsigned>()->default_value(0),
                    "Maximum number of iterations for the minimizer."  )
      ("maxpairs,m", po::value<LaDa::types::t_unsigned>()->default_value(5),
                     "Max distance for pairs (in neighbors)."  )
      ("minimizer", po::value<std::string>()->default_value("frprmn"),
                    "Type of Minimizer"  )
      ("iw", po::value<LaDa::types::t_real>()->default_value(0),
             "Initial value of the weights."  )
      ("loo", "Single leave-one-out procedure." )
      ("fit", "Single leave-none-out procedure." )
      ("alpha", po::value<LaDa::types::t_real>()->default_value(1),
                 "Lambda for pair regularization. With --loo/fit and --tcoef only. " )
      ("tcoef", po::value<LaDa::types::t_real>()->default_value(0),
                "\"t\" coefficient. With --loo/fit only. " )
      ("volkerreg", " Pair regulation done with Volker Blum's"
                    " normalization (changes the units of the \"t\" coefficient). " );
  leavemanyout.add_cmdl( specific );
     
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
 
  __BPO_MAP__
 
  std::cout << "\n" << __PROGNAME__
            << " from the " << PACKAGE_STRING << " package.\n";
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
  LADA_DO_NASSERT( not fs::exists( dir ), dir << " does not exist.\n" );
  LADA_DO_NASSERT( not fs::exists( dir / "LDAs.dat" ), 
              ( dir / "LDAs.dat" ) << " does not exist.\n" );
  LADA_DO_NASSERT
  (
    not (    fs::is_regular( dir / "LDAs.dat" ) 
          or fs::is_symlink( dir / "LDAs.dat" ) ),
    ( dir / "LDAs.dat" ) << " is neither a regular file nor a symlink.\n"
  ) 

  fs::path jtypes( vm["jtypes"].as< std::string >() );
  if(    ( not fs::exists( jtypes ) )
      or ( not ( fs::is_symlink(jtypes) or fs::is_regular(jtypes) ) ) )
    jtypes = dir / jtypes;
  LADA_DO_NASSERT( not fs::exists( jtypes ),
                 "Jtypes input file " << jtypes
              << " could not be found in ./ nor in " << dir << ".\n" )
  LADA_DO_NASSERT( not ( fs::is_regular( jtypes ) or fs::is_symlink( jtypes ) ),
              jtypes << " is a not a valid file.\n" );

  fs::path latinput( vm["latinput"].as< std::string >() );
  if(    ( not fs::exists( latinput ) )
      or ( not ( fs::is_symlink(latinput) or fs::is_regular(latinput) ) ) )
   latinput = dir / latinput;
  LADA_DO_NASSERT( not fs::exists( latinput ),
                 "Lattice input file " << latinput
              << " could not be found in ./ nor in " << dir << ".\n" )
  LADA_DO_NASSERT( not ( fs::is_regular( latinput ) or fs::is_symlink( latinput ) ),
              latinput << " is a not a valid file.\n" );

  const LaDa::types::t_unsigned maxpairs = vm["maxpairs"].as<LaDa::types::t_unsigned>();
  const LaDa::types::t_int verbosity = vm["verbose"].as<LaDa::types::t_int>();
  LaDa::types::t_unsigned seed = vm["seed"].as<LaDa::types::t_unsigned>();
  seed = LaDa::math::seed( seed );
  const LaDa::types::t_real tolerance( vm["tolerance"].as< LaDa::types::t_real >() );
  const LaDa::types::t_real linetol( vm["linetol"].as< LaDa::types::t_real >() );
  const LaDa::types::t_real linestep( vm["linestep"].as< LaDa::types::t_real >() );
  const std::string strategy( vm["strategy"].as< std::string >() );
  const LaDa::types::t_unsigned itermax( vm["itermax"].as< LaDa::types::t_unsigned >() );
  const std::string minimizer_type( vm["minimizer"].as< std::string >() );
  const LaDa::types::t_real iweights( vm["iw"].as< LaDa::types::t_real >() );
  const bool doloo = vm.count("loo") != 0;
  const bool dofit = vm.count("fit") != 0;
  const bool volkerreg = vm.count("volkerreg") != 0;
  const std::pair<LaDa::types::t_real, LaDa::types::t_real>
     pairreg( vm["alpha"].as<LaDa::types::t_real>(), vm["tcoef"].as<LaDa::types::t_real>() );

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
  // extract leave-many-out commandline
  leavemanyout.extract_cmdl( vm );
  std::cout << "\n";

  // consistency checks.
  LADA_NASSERT( LaDa::math::le(pairreg.second, 0e0), "Coefficient \"t\" cannot negative.\n" )

  // Loads lattice
  boost::shared_ptr< LaDa::Crystal::Lattice >
    lattice( LaDa::Crystal::read_lattice( latinput, dir ) );
  LaDa::Crystal::Structure::lattice = lattice.get();

  // Construct regularization.
  LaDa::CE::Regulated reg;
  reg.cgs.tolerance = 1e-5 * tolerance;
  reg.cgs.verbose = verbosity >= innermin;
  reg.cgs.itermax = 40;
  
  // reads jtypes
  LaDa::CE::Regulated :: t_Clusters clusters;
  // add pair terms.
  LaDa::CE::create_pairs( *lattice, maxpairs, clusters );
  std::cout << "Creating " << clusters.size() << " pair figures.\n";
  LaDa::CE::read_clusters( *lattice, jtypes, reg.clusters ); 
  std::cout << "Read " << reg.clusters.size()
            << " figures from " << jtypes << ".\n";
  reg.clusters.insert(   reg.clusters.begin()
                       + std::min( (size_t)2, reg.clusters.size() ), 
                       clusters.begin(), clusters.end() );
  clusters.clear();
  if( verbosity >= startclusters )
  {
    LaDa::CE::Regulated :: t_Clusters :: const_iterator i_class = reg.clusters.begin();
    LaDa::CE::Regulated :: t_Clusters :: const_iterator i_class_end = reg.clusters.end();
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
  if( doloo or dofit or leavemanyout.do_perform )
  {
    // Create object for fitting procedures.
    LaDa::CE::Fit< LaDa::CE::FittingPolicy::PairReg<> > fit;
    fit.alpha = pairreg.first;
    fit.tcoef = pairreg.second;
    fit.do_pairreg = ( not LaDa::math::is_null( pairreg.second ) and maxpairs );
    fit.laksreg = not volkerreg;
    fit.verbose = verbosity >= detailederrors;
    reg.cgs.verbose = verbosity >= outermin;
    // initialization
    LaDa::Crystal::read_ce_structures( dir / "LDAs.dat", fit.structures() );
    if( verbosity >= startstructures )
      std::for_each
      (
        fit.structures().begin(), fit.structures().end(), 
        std::cout << bl::_1 << "\n"
      );
    fit.init( reg.clusters );
    LaDa::opt::NErrorTuple nerror = fit.mean_n_var();
    std::cout << "Data mean: " << nerror.nmean() << "\n"
              << "Data variance: " << nerror.nvariance() << "\n"
              << "Data maximum error: " << nerror.nmax() << "\n";


    if( doloo )
    {
      LaDa::CE::BaseFit::t_Vector x( reg.clusters.size() );
      std::cout << "Starting Leave-One-Out Procedure.\n";
      std::pair< LaDa::opt::ErrorTuple, LaDa::opt::ErrorTuple > errors;
      errors = leave_one_out( fit, reg.cgs, x, verbosity >= serrors );
      std::cout << "Average Training Errors:\n " << ( nerror = errors.first ) << "\n";
      std::cout << "Final Prediction Errors:\n " << ( nerror = errors.second ) << "\n\n";
    }
    if( dofit )
    {
      LaDa::CE::BaseFit::t_Vector x( reg.clusters.size() );
      std::cout << "Starting Fitting Procedure.\n";
      LaDa::opt::ErrorTuple errors( fit( x, reg.cgs ) );
      std::cout << "Fitting Errors:\n " << ( nerror = errors ) << "\n\n";
      std::cout << "Fitting Errors (-log10):\n " << LaDa::opt::log( nerror = errors ) << "\n\n";
      if( verbosity >= printinteractionenergies )
      {
        LaDa::CE::BaseFit::t_Clusters :: const_iterator i_eclusters = reg.clusters.begin();
        LaDa::CE::BaseFit::t_Clusters :: const_iterator i_eclusters_end = reg.clusters.end();
        LaDa::CE::BaseFit::t_Vector :: const_iterator i_energ = x.begin();
        for(; i_eclusters != i_eclusters_end; ++i_eclusters, ++i_energ )
          std::cout << *i_energ << " " << i_eclusters->size()
                    <<  "\n" << i_eclusters->front() << "\n"; 
      }
    }
    if( leavemanyout.do_perform )
    {
      std::cout << "\nStarting leave-many out predictive fit." << std::endl;
      LaDa::opt::t_ErrorPair errors;
      leavemanyout.verbosity = verbosity - 1;
      errors = LaDa::CE::leave_many_out( leavemanyout, fit, reg.cgs );
      std::cout << "Average Training Errors:\n " << ( nerror = errors.first ) << "\n";
      std::cout << "Final Prediction Errors:\n "
                << ( nerror = errors.second ) << "\n\n";
    }
  }
  else
  {
    // initialization.
    LaDa::Crystal::read_ce_structures( dir / "LDAs.dat", reg.structures() );
    if( verbosity >= startstructures )
      std::for_each
      (
        reg.structures().begin(), reg.structures().end(), 
        std::cout << bl::_1 << "\n"
      );
    reg.init();
    LaDa::opt::NErrorTuple nerror = reg.mean_n_var();
  
    TiXmlElement fakexml( "Minimizer" );
    fakexml.SetAttribute( "type", minimizer_type );
    fakexml.SetDoubleAttribute( "tolerance", tolerance );
    fakexml.SetAttribute( "itermax", itermax );
    fakexml.SetDoubleAttribute( "linetolerance", linetol );
    fakexml.SetDoubleAttribute( "linestep", linestep );
    fakexml.SetAttribute( "strategy", strategy );
    fakexml.SetAttribute( "verbose", verbosity >= outermin ? "true": "false" );
#   if defined(LADA_WITH_GSL) or defined(LADA_WITH_MINUIT2)
    typedef LaDa::Minimizer::Frpr t_Minimizer;
#   else
    typedef LaDa::Minimizer::Variant
            < 
              boost::mpl::vector
              <
                LaDa::Minimizer::Frpr
#               ifdef LADA_WITH_GSL
                  , LaDa::Minimizer::Gsl
#               endif
#               ifdef LADA_WITH_MINUIT2
                  , LaDa::Minimizer::Minuit2
#               endif
              > 
            > t_Minimizer;
#   endif
    t_Minimizer minimizer;
    if( not minimizer.Load( fakexml ) )
    {
       std::cerr << "Could not load minimizer " << minimizer_type << ".\n";
       return 1;
    } 
    LaDa::CE::drautz_diaz_ortiz( reg, minimizer, verbosity, iweights);
  }
  return 0;
  __BPO_CATCH__()
}


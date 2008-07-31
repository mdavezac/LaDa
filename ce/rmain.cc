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
    po::options_description specific("Regulation Options");
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
      ("lambda", po::value<types::t_real>()->default_value(1),
                 "Lambda for pair regularization. With --loo only. " )
      ("tcoef", po::value<types::t_real>()->default_value(1),
                "\"t\" coefficient. With --loo only. " );
       
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
      std::cout << "Usage: " << argv[0] << " [options] DATADIR JTYPES\n"
                   "  _ DATATDIR (=./) is an optional path to the"
                   " training set input.\n"
                   "  _ JTYPES (=jtypes) is an optional filename"
                   " for the file \n"
                   "                 from which to load the lattice.\n"
                   "                 (DATADIR must be present when using JTYPES).\n"
                   "JTYPES should be\n"
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
    std::pair<types::t_real, types::t_real> pairreg( vm["lambda"].as<types::t_real>(),
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
    if( not doloo )
      std::cout << "   Non-linear Minimizer: " << minimizer_type << "\n"
                << "   Initial value of the weights " << iweights << "\n";
    else
      std::cout << "   lambda: " << pairreg.first << "\n"
                << "   t: " << pairreg.second << "\n";
    std::cout << "\n";

    // Loads lattice
    boost::shared_ptr< Crystal::Lattice >
      lattice( Crystal::read_lattice( latinput, dir ) );

    // read structures as lda\@nrel input.
    std::vector< Crystal :: Structure > structures;
    Crystal::read_ce_structures( dir / "LDAs.dat", structures );
    Crystal::Structure::lattice = lattice.get();
    if( verbosity >= 5 )
      std::for_each( structures.begin(),
                     structures.end(), 
                     std::cout << bl::_1 << "\n");

    // Construct regularization.
    CE::Regulated reg;
    reg.cgs.tolerance = tolerance;
    reg.cgs.verbose = verbosity >= 4;
    reg.cgs.itermax = 40;

    // reads jtypes
    CE::Regulated :: t_Clusters clusters;
    // add pair terms.
    CE::create_pairs( *lattice, maxpairs, clusters );
    std::cout << "Creating " << clusters.size()
              << " pair figures from " << jtypes << ".\n";
    CE::read_clusters( *lattice, jtypes, reg.clusters ); 
    CE::Regulated :: t_Clusters :: const_iterator i_class = clusters.begin();
    CE::Regulated :: t_Clusters :: const_iterator i_class_end = clusters.end();
    std::cout << "Read " << reg.clusters.size()
              << " figures from " << jtypes << ".\n";
    reg.clusters.insert( reg.clusters.begin() + std::min( (size_t)2, reg.clusters.size() ), 
                         clusters.begin(), clusters.end() );
    clusters.clear();
    if( verbosity >= 6 )
      for(; i_class != i_class_end; ++i_class )
      {
        std::for_each( i_class->begin(), i_class->end(), std::cout << bl::_1 << "\n" );
        std::cout << i_class->size() << "\n";
      }

    // initialization.
    reg.init( structures );
    


    // now for the real job.
    if( doloo ) // single leave-one-out with pair regulation.
    {
      CE::Regulated :: t_Vector weights( reg.clusters.size() );
      std::fill( weights.begin(), weights.end(), 0e0 );
      if( not Fuzzy::is_zero( pairreg.second ) and maxpairs )
      {
        CE::Regulated :: t_Clusters :: const_iterator i_clusters = reg.clusters.begin();
        CE::Regulated :: t_Clusters :: const_iterator i_clusters_end = reg.clusters.end();
        CE::Regulated :: t_Vector :: iterator i_weight = weights.begin();
         
        // Now loop over pairs.
        types::t_real normalization(0);
        for(; i_clusters != i_clusters_end; ++i_clusters, ++i_weight )
        {
          __ASSERT( i_clusters->size() > 0, "No clusters in class.\n" )
          if( i_clusters->front().size() != 2 ) continue;

          types::t_real D( i_clusters->size() );
          types::t_real R( atat::norm2(   i_clusters->front()[0] 
                                        - i_clusters->front()[1] ) );

          *i_weight = std::pow( R, pairreg.first ) / D;
          normalization += (*i_weight);
          *i_weight = std::sqrt( *i_weight );
        }
        normalization = std::sqrt( normalization * pairreg.second );
        std::for_each( weights.begin(), weights.end(),
                       bl::_1 *= bl::constant( normalization ) );
      }
      std::cout << "Leave-One-Out MSE: " << reg( &weights[0] ) << "\n"; 
    }
    else if( minimizer_type.compare( "simplex" ) == 0 )
    {
      Minimizer::Simplex simplex;
      simplex.tolerance = tolerance;
      simplex.verbose = verbosity >= 2;
      simplex.itermax = itermax;
      simplex.stepsize = 1;
      CE::drautz_diaz_ortiz( reg, simplex, verbosity - 1, iweights);
    }
    else 
    {
      Minimizer::Gsl gsl;
      gsl.type =  Minimizer::Gsl::SteepestDescent;
      gsl.tolerance = tolerance;
      gsl.verbose = verbosity >= 2;
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
    
      CE::drautz_diaz_ortiz( reg, gsl, verbosity - 1, iweights);
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

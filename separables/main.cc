#include "LaDaConfig.h"

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
#include <boost/algorithm/string.hpp>
#ifdef LADA_MPI
# include <boost/scoped_ptr.hpp>
# include <boost/mpi/environment.hpp>
#endif

#include <opt/mpi.h>
#include <minimizer/cgs.h>
#include <math/fuzzy.h>
#include <math/random.h>
#include <opt/types.h>
#include <opt/debug.h>
#include <opt/errors.h>
#include <opt/bpo_macros.h>
#include <crystal/lattice.h>
#include <crystal/structure.h>
#include <crystal/enumerate.h>
#include <ce/cluster.h>
#include <ce/regularization.h>
#include <ce/harmonic.h>

#include <opt/leave_many_out.h>

#include "functional.h"
#include "sepmappings.h"
#include "colmappings.h"
#include "collapse.h"
#include "prepare.h"
#include "methods.h"
#include "mixed.h"
#include "mixedfunctional.h"

#define __PROGNAME__ "Mixed Approach: CE and Sum of Separable functions on a fixed lattice." 
#if defined(_CUBIC_CE_)
typedef LaDa::CE::ConstituentStrain::Harmonic::Cubic t_Harmonic;
#elif defined( _TETRAGONAL_CE_ )
typedef LaDa::CE::ConstituentStrain::Harmonic::Tetragonal t_Harmonic;
#else
#error Please specify _CUBIC_CE_ or _TETRAGONAL_CE_
#endif

const LaDa::types::t_unsigned print_reruns   = 1;
const LaDa::types::t_unsigned print_checks   = 2;
const LaDa::types::t_unsigned print_function = 3;
const LaDa::types::t_unsigned print_allsq    = 4;
const LaDa::types::t_unsigned print_data     = 5;
const LaDa::types::t_unsigned print_llsq     = 6;

int main(int argc, char *argv[]) 
{
  namespace bl = boost::lambda;
  namespace fs = boost::filesystem;

  LADA_MPI_START
  LADA_TRY_BEGIN
  __BPO_START__
         ("verbose,p", po::value<LaDa::types::t_unsigned>()->default_value(0),
                       "Level of verbosity.\n"  )
         ("seed", po::value<LaDa::types::t_unsigned>()->default_value(0),
                  "Seed of the random number generator.\n"  );
  __BPO_SPECIFICS__( "Separables Options" )
      ("loo", "Performs leave-one-out cross-validation.\n"  )
      ("nofit", "Does not perform fit.\n"  )
      ("rank,r", po::value<LaDa::types::t_unsigned>()->default_value(3),
                 "Rank of the sum of separable functions." )
      ("basis",  po::value<std::string>()->default_value("1x1x4"),
                 "Description of the ranks/size of the figure used\n." )
      ("tolerance", po::value<LaDa::types::t_real>()->default_value(1e-4),
                    "Tolerance of the alternating linear-least square fit.\n"  )
      ("maxiter", po::value<LaDa::types::t_unsigned>()->default_value(40),
                  "Maximum number of iterations for"
                  " Alternating linear-least square fit.\n"  )
      ("1dtolerance", po::value<LaDa::types::t_real>()->default_value(1e-4),
                      "Tolerance of the 1d linear-least square fit.\n" ) 
      ("random", po::value<LaDa::types::t_real>()->default_value(5e-1),
                 "Coefficients' randomness.\n" )
      ("lambda,l", po::value<LaDa::types::t_real>()->default_value(0),
                   "Regularization factor.\n" )
      ("maxpairs,m", po::value<LaDa::types::t_unsigned>()->default_value(5),
                     "Max distance for pairs (in neighbors)."  )
      ("J0", "Introduce J0 in clusters."  )
      ("J1", "Introduce J1 in clusters."  )
      ("rm", "Remove clusters which are contained within the separable-basis."  )
      ("alpha", po::value<LaDa::types::t_real>()->default_value(1),
                 "Lambda for pair regularization. With --loo/fit and --tcoef only. " )
      ("tcoef", po::value<LaDa::types::t_real>()->default_value(0),
                "\"t\" coefficient. With --loo/fit only. " )
      ("volkerreg", " Pair regulation done with Volker Blum's"
                    " normalization (changes the units of the \"t\" coefficient). " )
      ("bestof,b", po::value<LaDa::types::t_unsigned>()->default_value(1),
                   "Performs best-of fit.\n" )
      ("enum", po::value<std::string>(), "Enumerate PI-file.\n" )
      ("genes", po::value<std::string>()->default_value(""), "Figure bitstring.\n" )
      ("jtypes", po::value<std::string>(), "Figure description file.\n" )
      ("cs", po::value<std::string>(), "Constituent Strain input file.\n" )
      ("which,w", po::value<LaDa::types::t_unsigned>()->default_value(0),
                   "Performs best-of for 0 (variance), 1 (mean), or 2(max).\n" )
      ("gus", "Use Gus Harts' type PI-file for structural enumeration." )
      ("print", po::value<std::string>()->default_value(""),
                   "Prints out: \"function\".\n" );
  LaDa::Fitting::LeaveManyOut leavemanyout;
  leavemanyout.add_cmdl( specific );
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
 
  std::cout << "\n" << __PROGNAME__ \
            << " from the " << PACKAGE_STRING << " package.\n";
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
  LADA_DO_NASSERT( not fs::exists( dir ), dir << " does not exist.\n" );
  LADA_DO_NASSERT( not fs::exists( dir / "LDAs.dat" ), 
              ( dir / "LDAs.dat" ) << " does not exist.\n" );
  LADA_DO_NASSERT
  (
    not (    fs::is_regular( dir / "LDAs.dat" ) 
          or fs::is_symlink( dir / "LDAs.dat" ) ),
    ( dir / "LDAs.dat" ) << " is neither a regular file nor a symlink.\n"
  ) 
  fs::path latinput( vm["latinput"].as< std::string >() );
  if(    ( not fs::exists( latinput ) )
      or ( not ( fs::is_symlink(latinput) or fs::is_regular(latinput) ) ) )
   latinput = dir / latinput;
  LADA_DO_NASSERT( not fs::exists( latinput ),
                 "Lattice input file " << latinput
              << " could not be found in ./ nor in " << dir << ".\n" )
  LADA_DO_NASSERT( not ( fs::is_regular( latinput ) or fs::is_symlink( latinput ) ),
              latinput << " is a not a valid file.\n" );

  fs::path doenum("");
  if( vm.count( "enum" ) > 0  )
  {
    doenum = vm["enum"].as< std::string >();
    fs::path orig = doenum;
    if(    ( not fs::exists( doenum ) )
        or ( not ( fs::is_symlink(doenum) or fs::is_regular(doenum) ) ) )
      doenum = dir / doenum;
    LADA_DO_NASSERT( not fs::exists( doenum ),
                   "Lattice input file " << orig
                << " could not be found in ./ nor in " << doenum << ".\n" )
    LADA_DO_NASSERT( not ( fs::is_regular( doenum ) or fs::is_symlink( doenum ) ),
                doenum << " is a not a valid file.\n" );
  }
  fs::path jtypes("");
  if( vm.count( "jtypes" ) > 0  )
  {
    jtypes = vm["jtypes"].as< std::string >();
    fs::path orig = jtypes;
    if(    ( not fs::exists( jtypes ) )
        or ( not ( fs::is_symlink(jtypes) or fs::is_regular(jtypes) ) ) )
      jtypes = dir / jtypes;
    LADA_DO_NASSERT( not fs::exists( jtypes ),
                   "Figure descrition file " << orig
                << " could not be found in ./ nor in " << jtypes << ".\n" )
    LADA_DO_NASSERT( not ( fs::is_regular( jtypes ) or fs::is_symlink( jtypes ) ),
                jtypes << " is a not a valid file.\n" );
  }
  fs::path cs("");
  if( vm.count( "cs" ) > 0 and vm.count( "enum" ) > 0  )
  {
    cs = vm["cs"].as< std::string >();
    fs::path orig = cs;
    if(    ( not fs::exists( cs ) )
        or ( not ( fs::is_symlink(cs) or fs::is_regular(cs) ) ) )
      cs = dir / cs;
    LADA_DO_NASSERT( not fs::exists( cs ),
                   "Constituent-strain input file " << orig
                << " could not be found in ./ nor in " << cs << ".\n" )
    LADA_DO_NASSERT( not ( fs::is_regular( cs ) or fs::is_symlink( cs ) ),
                cs << " is a not a valid file.\n" );
  }

  const bool dofit( vm.count("nofit") == 0 );
  const bool doloo( vm.count("loo") != 0 );
  const LaDa::types::t_unsigned verbosity = vm["verbose"].as<LaDa::types::t_unsigned>();
  LaDa::types::t_unsigned seed = vm["seed"].as<LaDa::types::t_unsigned>();
  seed = LaDa::math::seed( seed );
  const LaDa::types::t_unsigned rank( vm["rank"].as< LaDa::types::t_unsigned >() );
  const std::string bdesc( vm["basis"].as<std::string>() );
  const LaDa::types::t_real tolerance( vm["tolerance"].as< LaDa::types::t_real >() );
  const LaDa::types::t_unsigned maxiter( vm["maxiter"].as< LaDa::types::t_unsigned >() );
  const LaDa::types::t_real dtolerance( vm["1dtolerance"].as< LaDa::types::t_real >() );
  const LaDa::types::t_real howrandom( vm["random"].as<LaDa::types::t_real>() );
  const LaDa::types::t_real lambda( vm["lambda"].as<LaDa::types::t_real>() );
  const bool volkerreg = vm.count("volkerreg") != 0;
  const LaDa::types::t_real alpha( vm["alpha"].as<LaDa::types::t_real>() );
  const LaDa::types::t_real tcoef( vm["tcoef"].as<LaDa::types::t_real>() );
  const LaDa::types::t_unsigned maxpairs = vm["maxpairs"].as<LaDa::types::t_unsigned>();
  const bool J0( vm.count("J0") > 0 );
  const bool J1( vm.count("J1") > 0 );
  const bool gusfile( vm.count("gus") > 0 );
  const bool rmpairs( vm.count("rm") > 0 );
  LADA_NASSERT( LaDa::math::le(tcoef, 0e0), "Coefficient \"t\" cannot negative.\n" )
  const LaDa::types::t_unsigned bestof( vm["bestof"].as<LaDa::types::t_unsigned>() );
  LADA_DO_NASSERT( bestof == 0, "0 jobs to be performed..." )
  const LaDa::types::t_unsigned which( vm["which"].as<LaDa::types::t_unsigned>() );
  const std::string genes( boost::algorithm::trim_copy( vm["genes"].as<std::string>() ) );
  if( not genes.empty() )
    LADA_DO_NASSERT( not boost::algorithm::all( genes, boost::algorithm::is_any_of( "01" ) ),
                "Unknown bitstring format.\n" )
  LADA_DO_NASSERT( which >= 3, "Don't know which error to perform bestof for.\n" )
  const std::string print( vm["print"].as<std::string>() );

  if( not ( doenum.empty() or dofit ) )
    std::cout << "Must perform a fit when attempting PIfile enumeration.\n";

  // Loads lattice
  boost::shared_ptr< LaDa::Crystal::Lattice >
    lattice( LaDa::Crystal::read_lattice( latinput, dir ) );
  LaDa::Crystal::Structure::lattice = lattice.get();

  // create pair terms.
  typedef std::vector< std::vector< LaDa::CE::Cluster > > t_Clusters;
  std::vector< std::vector< LaDa::CE::Cluster > > clusters;
  LaDa::CE :: create_pairs( *lattice, maxpairs, clusters );
  std::cout << "  Creating " << clusters.size() << " pair figures.\n";
  if( not jtypes.empty() )
  {
    const size_t d( clusters.size() );
    LaDa::CE::read_clusters( *lattice, jtypes, clusters, genes ); 
    std::cout << "  Read " << clusters.size() - d
              << " from input file " << jtypes << "\n";
  }
  LaDa::CE::Cluster cluster;
  if( J0 )
  {
    bool isfound = false;
    foreach( t_Clusters :: value_type &_class, clusters )
      if( _class.front().size() == 0 ) { isfound = true; break; }
    if( isfound ) clusters.push_back( std::vector<LaDa::CE::Cluster>(1, cluster) ); 
  }
  cluster.Vectors().resize(1, LaDa::math::rVector3d(0,0,0) );
  if( J1 )
  {
    bool isfound = false;
    foreach( t_Clusters :: value_type &_class, clusters )
      if( _class.front().size() == 1 ) { isfound = true; break; }
    if( isfound ) clusters.push_back( std::vector<LaDa::CE::Cluster>(1, cluster) ); 
  }
  
  // Initializes fitting.
  typedef LaDa::Fitting::AlternatingLeastSquare<LaDa::Fitting::Cgs> t_Fitting;
  t_Fitting allsq;
  allsq.itermax = maxiter;
  allsq.tolerance = tolerance;
  allsq.verbose = verbosity >= print_allsq;
  allsq.linear_solver.tolerance = dtolerance;
  allsq.linear_solver.verbose = verbosity >= print_llsq;

  // Initializes basis.
  LaDa::CE::PosToConfs postoconfs( *LaDa::Crystal::Structure::lattice );
  if(  ( not bdesc.empty() ) and rank )
    postoconfs.create_positions( bdesc );

  // Separables traits.
  typedef LaDa::Traits::CE::Separables< LaDa::CE::Mapping::VectorPlus<2> > t_FunctionTraits;
  typedef LaDa::CE::Separables< t_FunctionTraits > t_Function;
  // Collapse Traits
  typedef LaDa::CE::Mapping::SymEquiv t_Mapping;
  typedef LaDa::CE::Policy::Regularization< t_Function > t_Regularization;
  typedef boost::numeric::ublas::matrix<size_t> t_Confs;
  typedef LaDa::CE::Policy::HighMemUpdate< t_Function, t_Mapping, t_Confs > t_UpdatePolicy;
  typedef LaDa::Traits::CE::Collapse< t_Function, t_Mapping, 
                                      t_Regularization, t_Confs,
                                      t_UpdatePolicy > t_CollapseTraits;
  // CE base
  typedef LaDa::CE::Fit< LaDa::CE::FittingPolicy::PairReg<> > t_CEBase;
  // Mixed approach traits.
  typedef LaDa::Traits::CE::MixedApproach< LaDa::CE::Collapse<t_CollapseTraits>, 
                                           t_CEBase > t_MixedTraits;
  typedef LaDa::CE::MixedApproach< t_MixedTraits > t_Collapse;

  // Finally creates mixed approach object.
  t_Collapse mixed;
  // initializes ce part.
  LaDa::Crystal::read_ce_structures( dir / "LDAs.dat", mixed.cefit().structures() );
  if( verbosity >= print_data )
    std::for_each
    (
      mixed.cefit().structures().begin(), mixed.cefit().structures().end(), 
      std::cout << bl::_1 << "\n"
    );
  mixed.cefit().alpha = alpha;
  mixed.cefit().tcoef = tcoef;
  mixed.cefit().do_pairreg = ( not LaDa::math::is_null( alpha ) and maxpairs );
  mixed.cefit().laksreg = not volkerreg;
  mixed.cefit().verbose = verbosity >= print_checks;
  // initializes collapse part.
  mixed.collapse().init( mixed.cefit().structures(), postoconfs );
  mixed.collapse().regularization().lambda = lambda;
  // initializes mixed.
  std::copy( clusters.begin(), clusters.end(),
             std::back_inserter( mixed.clusters() ) );
  if( verbosity >= print_data )
  {
    LaDa::CE::Regulated :: t_Clusters :: const_iterator i_class = mixed.clusters().begin();
    LaDa::CE::Regulated :: t_Clusters :: const_iterator i_class_end =
       mixed.clusters().end();
    for(; i_class != i_class_end; ++i_class )
    {
      if( verbosity < print_data + 1  )  std::cout << i_class->front()
                                                   << " D=" << i_class->size() 
                                                   << "\n";
      else
      {
        std::cout << "Cluster Class: " << i_class->size() << "\n";
        std::for_each( i_class->begin(), i_class->end(), 
                       std::cout << bl::_1 << "\n" );
      }
    }
  }

  if( rmpairs and postoconfs.positions.size() )
  {
    LaDa::CE::remove_contained_clusters( postoconfs.positions, mixed.clusters() );
    size_t size( mixed.clusters().size() );
    if( J0 ) --size;
    if( J1 ) --size;
    std::cout << "Retained " << size << " pair figures.\n";
  }
  clusters.clear();
  mixed.init( rank, postoconfs.dof() );

  LaDa::opt::NErrorTuple nerror( LaDa::opt::mean_n_var( mixed.cefit().structures() ) ); 

  std::cout << "Shape of separable function: " << bdesc << "\n"
            << "Rank of separable function " << rank << "\n"
            << "Size of separable function "
            << postoconfs.dof() << "\n"
            << "Data directory: " << dir << "\n"
            << "Number of data points: " << mixed.mapping().size() << " for " 
            << mixed.nbconfs()
            << " symetrically inequivalent configurations.\n";
  if( not verbosity ) std::cout << "Quiet output.\n";
  else std::cout << "Level of verbosity: " << verbosity << "\n";
  std::cout << "Alternating linear-least square tolerance: " 
               << tolerance << "\n"
            << "Maximum number of iterations for alternating least-square fit: "
               << maxiter << "\n"
            << "1d linear-least square tolerance: " 
               << dtolerance << "\n"
            << "Data mean: " << nerror.nmean() << "\n"
            << "Data variance: " << nerror.nvariance() << "\n"
            << "Data maximum error: " << nerror.nmax() << "\n"
            << "Random Seed: " << seed << "\n"
            << "Separables Regulation factor: " << lambda << "\n"
            << "Maximum distance of pairs (nth neighbor): " << maxpairs << "\n";
  if( rmpairs ) std::cout << "Will remove clusters completely covered"
                             " by the separable-function.\n";
  if( J0 ) std::cout << "Introducing J0 in figures.\n";
  if( J1 ) std::cout << "Introducing J1 in figures.\n";
  std::cout << "CE regulation alpha: " << alpha << "\n"
            << "CE t-coef: " << tcoef << "\n";
  if( volkerreg ) std::cout << "   Using Volker's normalization for t.\n";
  std::cout << "Performing best of  " << bestof;
  if( which == 0 ) std::cout << " for variance\n";
  if( which == 1 ) std::cout << " for mean\n";
  if( which == 2 ) std::cout << " for max\n";
  if( not ( cs.empty() and doenum.empty() ) )
    std::cout << "Constituent-strain input file: " << cs << "\n";
  if( not jtypes.empty() ) std::cout << "Figures description file: " << jtypes << "\n";
  if( not doenum.empty() ) std::cout << "Enumeration PI-file: " << doenum << "\n";
  if( not genes.empty() ) std::cout << "genes: " << genes << "\n";
  // extract leave-many-out commandline
  leavemanyout.extract_cmdl( vm );

  // Initializes best of fit.
  LaDa::CE::Method::Fit< LaDa::CE::Method::Policy::BestOf< LaDa::CE :: CollapseState > >
    fit( mixed.cefit().structures() );
  fit.verbosity = verbosity -1;
  fit.policy.restarts = bestof;
  fit.policy.which = which;
  fit.policy.howrandom = howrandom;

  // fitting.
  if( doloo )
  {
    typedef LaDa::CE::Mapping::ExcludeOne< t_Collapse::t_Traits::t_Mapping > t_looMapping;
    typedef LaDa::Traits::CE::MixedApproachWithNewMapping
            < 
              t_Collapse,
              t_looMapping 
            > :: type t_looCollapse;
    t_looCollapse loomixed; loomixed.operator=( mixed );
    loomixed.collapse().mapping().do_exclude = true;
    std::cout << "Starting Leave-One-Out Procedure.\n";
    LaDa::opt::t_ErrorPair errors;
    errors = LaDa::CE::Method::leave_one_out( loomixed, fit, allsq, verbosity - 1 );
    std::cout << "Average Training Errors:\n " << ( nerror = errors.first ) << "\n";
    std::cout << "Final Prediction Errors:\n " << ( nerror = errors.second ) << "\n\n";
  }
  if( leavemanyout.do_perform )
  {
    typedef std::vector< LaDa::types::t_unsigned > t_Excluded;
    typedef LaDa::CE::Mapping::ExcludeMany
            < 
              t_Collapse :: t_Traits :: t_Mapping, 
              t_Excluded
            > t_lmoMapping;
    typedef LaDa::Traits::CE::MixedApproachWithNewMapping
            < 
              t_Collapse,
              t_lmoMapping 
            > :: type t_lmoCollapse;
    t_lmoCollapse lmomixed; lmomixed.operator=( mixed );
    std::cout << "\nStarting leave-many out predictive fit." << std::endl;
    LaDa::Fitting::LeaveManyOut::t_Return errors;
    leavemanyout.verbosity = verbosity - 1;
    errors = LaDa::CE::Method::leave_many_out( leavemanyout, lmomixed, fit, allsq );
    std::cout << "Average Training Errors:\n " << ( nerror = errors.first ) << "\n";
    std::cout << "Final Prediction Errors:\n "
              << ( nerror = errors.second ) << "\n\n";
    if( print.find("function") != std::string::npos )
    {
      lmomixed.reassign();
      std::cout << lmomixed << "\n";
    }
  }
  if( dofit or (not doenum.empty()) )
  {
    std::cout << "\nFitting using whole training set:" << std::endl;
    fit.verbosity = verbosity;
    nerror = fit( mixed, allsq );
    mixed.reassign();
    std::cout << nerror << "\n"; 
    if( verbosity >= print_function or print.find("function") != std::string::npos )
      std::cout << mixed << "\n";
  }
  if( not doenum.empty() )
  {
    std::cout << "\nStarting Exhaustive Enumeration of " << doenum << "\n\n";
    typedef LaDa::CE::MixedSeparables
            < 
              t_Collapse :: t_Traits :: t_Separables, 
              t_Harmonic 
            > t_MixedFunctional;
    t_MixedFunctional mixedfunc( *lattice );
    mixedfunc.init( cs.string(), bdesc );
    mixedfunc = mixed;
    LaDa::Crystal :: Structure structure;
    structure.cell = lattice->cell;
    size_t n(0);
    foreach( const LaDa::Crystal::Lattice::t_Site &site, lattice->sites )
    {
      LaDa::Crystal::Structure::t_Atom atom;
      atom.pos = site.pos;
      atom.site = n; ++n;
      atom.type = -1e0;
      structure.atoms.push_back( atom );
    }
    structure.find_k_vectors();
    std::cout << "  @-0 " << structure.get_concentration()
              << " " << mixedfunc( structure ) << "\n";
    if( not gusfile ) LaDa::Crystal::enumerate_pifile( doenum.string(), mixedfunc );
    else LaDa::Crystal::enumerate_gusfile( doenum.string(), mixedfunc );
    foreach( LaDa::Crystal::Structure::t_Atom &atom, structure.atoms )
      atom.type *= -1e0;
    std::cout << "  @0 " << structure.get_concentration()
              << " " << mixedfunc( structure ) << "\n";
  }
  std::cout << "\n\n\nEnd of " << __PROGNAME__ << ".\n" << std::endl;

  LADA_MPI_CODE( MPI_Finalize() );
  return 1;
  __BPO_CATCH__( LADA_MPI_CODE( MPI_Finalize() ) )
}


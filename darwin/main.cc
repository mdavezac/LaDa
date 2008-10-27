//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdexcept>       // std::runtime_error
#ifdef _MPI
# include <boost/mpi/environment.hpp>
#endif
#include <boost/scoped_ptr.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <opt/initial_path.h>



#ifdef _PESCAN
# include "bandgap.h"
  typedef BandGap :: Evaluator t_Evaluator;
# define __PROGNAME__ "Band-Gap Optimization"
#elif defined(_CE)
# include "groundstate.h"
  typedef GroundState :: Evaluator t_Evaluator;
# define __PROGNAME__ "Cluster Expansion Optimization"
#elif defined(_MOLECULARITY)
# include "molecularity.h"
  typedef Molecularity :: Evaluator t_Evaluator;
# define __PROGNAME__ "Band-Gap Optimization for Epitaxial Structure"
#elif defined(_EMASS)
# include "emass.h"
  typedef eMassSL :: Evaluator t_Evaluator;
# define __PROGNAME__ "emass_opt"
#elif defined( _ALLOY_LAYERS_ )
# include "two_sites.h"
# include "alloylayers/evaluator.h"
# include "alloylayers/policies.h"
# include "alloylayers/object.h"
  typedef Individual :: Types
          < 
            GA :: AlloyLayers :: Object,
            TwoSites :: Concentration,
            TwoSites :: Fourier 
          > :: t_Vector t_Individual;
  typedef GA :: AlloyLayers :: Evaluator
          <
            t_Individual,
            GA :: AlloyLayers :: Translate,
            GA :: AlloyLayers :: AssignSignal
          > t_Evaluator;
# define __PROGNAME__ "alloylayers_opt"
#else 
# error Need to define _CE or _PESCAN or _MOLECULARITY or _ALLOY_LAYERS_
#endif

#include <revision.h>
#include <print/xmg.h>
#include <print/stdout.h>
#include <opt/debug.h>
#include <mpi/mpi_object.h>
#include <opt/bpo_macros.h>
#include <opt/initial_path.h>

#include "individual.h"
#include "darwin.h"


int main(int argc, char *argv[]) 
{
  namespace fs = boost::filesystem;
  namespace bl = boost::lambda;
  opt::InitialPath::init();

  __MPI_START__
  __TRYBEGIN

  __BPO_START__;
  __BPO_SPECIFICS__( "GA Options" )
    __BPO_RERUNS__;
  __BPO_GENERATE__()
  __BPO_MAP__
  __BPO_HELP_N_VERSION__


  fs::path input( vm["input"].as< std::string >() );
  __DOASSERT( not ( fs::is_regular( input ) or fs::is_symlink( input ) ),
              input << " is a not a valid file.\n" );
  const unsigned reruns = vm["reruns"].as< unsigned >();
  
  __ROOTCODE
  (
    (*::mpi::main),
    std::cout << "Will load input from file: " << input << ".\n"
              << "Will perform " << reruns << " GA runs.\n\n";
  )
    
  __MPICODE
  (
    std::string input_str = input.string();
    boost::mpi::broadcast( (*::mpi::main), input_str, 0 ); 
    input = input_str;
  )
  
  
  for( types::t_int i = 0; i < reruns; ++i )
  {
    GA::Darwin< t_Evaluator > ga;
    Print :: out << "load result: " << ga.Load(input.string()) << Print::endl; 
   //__DOASSERT( not ga.Load(input.string()),
   //            "Could not load input from file " << input << "\n" )
#   ifdef _ALLOY_LAYERS_
      ga.evaluator.do_dipole = true;
      ga.evaluator.connect
      (
        bl::_2[0] = bl::bind(  &GA::AlloyLayers::inplane_stress< t_Evaluator >,
                               bl::_1, bl::constant(ga.evaluator) ) 
      );
      ga.evaluator.connect
      (
        bl::_2[1] = bl::bind(  &GA::Keepers::OscStrength::osc_strength, bl::_1 )
      );
#   endif
    
    Print::out << "Rerun " << i+1 << " of " << reruns << Print::endl;
    __ROOTCODE
    (
      (*::mpi::main),
      std::cout << "Rerun " << i+1 << " of " << reruns << ".\n";
    )
    ga.run();
    // for reruns.
    Print::xmg.dont_truncate();
    Print::out.dont_truncate();
  }
  
  __BPO_CATCH__

  return 0;
}

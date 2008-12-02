//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdexcept>       // std::runtime_error
#include <functional>
#ifdef _MPI
# include <boost/mpi/environment.hpp>
#endif
#include <boost/scoped_ptr.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/lambda/bind.hpp>



// Includes optimization specific stuff.
# include "groundstates/main.extras.h"
# include "alloylayers/main.extras.h"

#include <revision.h>
#include <print/xmg.h>
#include <print/stdout.h>
#include <opt/debug.h>
#include <mpi/mpi_object.h>
#include <opt/bpo_macros.h>
#include <opt/initial_path.h>
#include <factory/factory.h>
#include <checkpoints/checkpoint.h>
#include <checkpoints/print_iteration.h>

#include "individual.h"
#include "darwin.h"
#include "taboos/populations.h"
#include "operators/populate.h"
#include "checkpoints/maxgen.h"
#include "checkpoints/maxeval.h"
#include "checkpoints/print_populations.h"
#include "checkpoints/stop_onfile.h"
#include "checkpoints/average_fitness.h"
#include "checkpoints/true_census.h"
#include "checkpoints/xcrysdenanim.h"
#include "checkpoints/xyzanim.h"



int main(int argc, char *argv[]) 
{
  namespace fs = boost::filesystem;
  namespace bl = boost::lambda;
  LaDa::opt::InitialPath::init();

  __MPI_START__
  __TRYBEGIN

  __BPO_START__;
  __BPO_HIDDEN__;
  __BPO_SPECIFICS__( "GA Options" )
    // optimization extra parameters.
#   include "groundstates/main.extras.h" 
#   include "alloylayers/main.extras.h" 
    __BPO_RERUNS__;
  __BPO_GENERATE__()
  __BPO_MAP__

  __BPO_PROGNAME__
  __BPO_VERSION__

  fs::path input( vm["input"].as< std::string >() );
  const unsigned reruns = vm["reruns"].as< unsigned >();
  std::string short_description = "";
  

  // Reads program options for alloylayers.
  // Prints program options to standard output.
#   include "alloylayers/main.extras.h" 
#   include "groundstates/main.extras.h" 
  if ( vm.count("help") ) 
  { 
    LaDa::GA::Darwin< t_Evaluator > ga;
#   include "main.extras.h" 
#   include "alloylayers/main.extras.h" 
#   include "groundstates/main.extras.h" 
    __ROOTCODE( (*::LaDa::mpi::main), 
      std::cout << "Description: " << short_description
                << "Usage: " << argv[0] << " [options] file.xml\n" 
                << "  file.xml is an optional filename for XML input.\n" 
                << "  Default input is input.xml.\n\n" 
                << all << "\n"
                << "The GA tag accepts the following attributes:\n"
                << ga.att_factory << "\n"
                << "The following GA Operators can be used in <Breeding>...</Breeding>:\n"
                << ga.operator_factory
                << "\nThe taboos can be chosen from the amongst the following:\n"
                << ga.taboo_factory
#     ifdef _ALLOY_LAYERS_
                << "\nThe following physical properties can be chosen for optimization:\n"
                  << properties_factory 
                  << "when included as in\n"
                     "  <Objectives>\n"
                     "    <Objective value=\"property A\" ... >\n"
                     "    <Objective value=\"property B\" ... >\n"
                     "    ....\n"
                     "  </Objectives>\n"
                     "Mind the plural and singulars ;).\n\n"
#     else          
                << "\n"
#     endif
                << "Miscellaeneous Items:"
                << ga.checkpoint_factory << "\n";
    ) 
    return 1; 
  } 

  __DOASSERT( not ( fs::is_regular( input ) or fs::is_symlink( input ) ),
              input << " is a not a valid file.\n" );
 

  __ROOTCODE
  (
    (*::LaDa::mpi::main),
    std::cout << "Will load input from file: " << input << ".\n";
  )
    
  __ROOTCODE
  (
    (*::LaDa::mpi::main),
    std::cout << "Will perform " << reruns << " GA runs.\n\n";
  )
  
  for( LaDa::types::t_int i = 0; i < reruns; ++i )
  {
    LaDa::GA::Darwin< t_Evaluator > ga;
    // creates factories
#   include "main.extras.h" 
#   include "groundstates/main.extras.h" 
#   include "alloylayers/main.extras.h" 
    // Connects assignement and print functors.
#   include "alloylayers/main.extras.h" 
#   include "groundstates/main.extras.h" 
    // Then loads input.
    LaDa::Print :: out << "load result: "
                       << ga.Load(input.string()) 
                       << LaDa::Print::endl; 

    
    LaDa::Print::out << "Rerun " << i+1 << " of " << reruns << LaDa::Print::endl;
    __ROOTCODE
    (
      (*::LaDa::mpi::main),
      std::cout << "Rerun " << i+1 << " of " << reruns << ".\n";
    )
    ga.run();
    // for reruns.
    LaDa::Print::xmg.dont_truncate();
    LaDa::Print::out.dont_truncate();
  }
  
  __BPO_CATCH__( __MPICODE( MPI_Finalize() ) )

  return 0;
}

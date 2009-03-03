//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdexcept>       // std::runtime_error
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

#define __PROGNAME__ "coollj"


#include <revision.h>
#include <opt/debug.h>
#include <opt/bpo_macros.h>
#include <opt/tuple_io.h>
#include <mpi/macros.h>
#ifdef _MPI
# include <boost/mpi/environment.hpp>
#endif

#include <crystal/structure.h>
#include <crystal/read_poscar.h>
#include <crystal/fractional_cartesian.h>
#include "clj.h"


int main(int argc, char *argv[]) 
{
  namespace fs = boost::filesystem;
  namespace bl = boost::lambda;

  __MPI_START__
  __TRYBEGIN

  __BPO_START__;
  po::options_description hidden( "hidden" ); 
  hidden.add_options()
    ("input, i", po::value<std::string>()->default_value("ep.input"),
     "Functional input filename" )
    ("poscar, p", po::value<std::string>()->default_value("POSCAR_0"),
     "POSCAR input filename" );
  __BPO_SPECIFICS__( "Coulomb + Lennard-Jones functional" );
  __BPO_GENERATE__()
  p.add("poscar", 2); 
  __BPO_MAP__
  __BPO_PROGNAME__
  __BPO_VERSION__
  if( vm.count("help") )
  {
    __MPICODE( boost::mpi::communicator world; )
    __ROOTCODE( world, \
      std::cout << argv[0] << " runs the Coulom+Lennard-Jones functional.\n" 
                   "Usage: " << argv[0] << " [options] FunctionalInput POSCAR\n" 
                   "  Defaults are ep.input and POSCAR_0, respectively.\n"
                << all << "\n";
    ) 
    return 0; 
  }


  // filenames.
  fs::path input( vm["input"].as< std::string >() );
  __DOASSERT( not ( fs::is_regular( input ) or fs::is_symlink( input ) ),
              input << " is a not a valid file.\n" );
  fs::path poscar( vm["poscar"].as< std::string >() );
  __DOASSERT( not ( fs::is_regular( poscar ) or fs::is_symlink( poscar ) ),
              poscar << " is a not a valid file.\n" );
  
  // data declaration.
  std::vector<std::string> species;
  LaDa::Models::Clj clj;
  LaDa::Crystal::TStructure<std::string> structure, forces;

  // loads stuff.
  LaDa::Models::read_fortran_input( clj, species, input );
  LaDa::Crystal::read_poscar( structure, poscar, species );
  forces = structure;
  
  // action.
  LaDa::Crystal::to_fractional( structure );
  const LaDa::types::t_real energy( clj( structure, forces ) );
// LaDa::Crystal::to_cartesian( structure );

  // printout.
  structure.energy = energy;
  std::cout << "Energy: " << energy  << "\n"
            << structure << "\n" << forces << "\n";
  __BPO_CATCH__( __MPICODE( MPI_Finalize() ) )

  return 0;
}

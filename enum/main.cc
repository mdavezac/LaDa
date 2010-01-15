//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdexcept>       // std::runtime_error
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

#define __PROGNAME__ "enum"


#include <opt/debug.h>
#include <opt/bpo_macros.h>
#include <opt/tuple_io.h>
#include <opt/random.h>
#include <mpi/mpi_object.h>
#ifdef _MPI
# include <boost/mpi/environment.hpp>
#endif

#include <crystal/structure.h>
#include <crystal/lattice.h>
#include "find_all_cells.h"


namespace Crystal = LaDa :: Crystal;

int main(int argc, char *argv[]) 
{
  namespace fs = boost::filesystem;
  namespace bl = boost::lambda;

  __MPI_START__
  __TRYBEGIN

  __BPO_START__;
  __BPO_HIDDEN__;
  __BPO_SPECIFICS__( "Enumeration" )
    // extra parameters for alloy layers.
    ("n", po::value<size_t>()->default_value(3), "Max size.\n" );
  __BPO_GENERATE__()
  __BPO_MAP__
  __BPO_PROGNAME__
  __BPO_VERSION__
  if( vm.count("help") )
  {
    __MPICODE( boost::mpi::communicator world; )
    __ROOTCODE( world, \
      std::cout << argv[0] << " Enumerates structures.\n" 
                << all << "\n";
    ) 
    return 0; 
  }

  fs::path input( vm["input"].as< std::string >() );
  __DOASSERT( not ( fs::is_regular( input ) or fs::is_symlink( input ) ),
              input << " is a not a valid file.\n" );
  
  boost::shared_ptr< LaDa::Crystal::Lattice >
    lattice( LaDa::Crystal::read_lattice( input ) );
  LaDa::Crystal::Structure::lattice = lattice.get();
  size_t N( vm["n"].as<size_t>() );

  boost::shared_ptr< std::vector<LaDa::Eigen::Matrix3d> > cells
    = LaDa::enumeration::find_all_cells(*lattice, N);

  __BPO_CATCH__( __MPICODE( MPI_Finalize() ) )

  return 0;
}

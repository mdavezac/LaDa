#include "LaDaConfig.h"

#include <boost/program_options.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>
#include <stdexcept>

#include <tinyxml/tinyxml.h>

#include <crystal/lattice.h>
#include <crystal/structure.h>
#include <crystal/fourier.h>
#include <opt/debug.h>
#include <opt/bpo_macros.h>
#include <opt/tuple_io.h>

#ifdef _LAYERED
# include <vff/va.h>
# include <vff/layered.h>
  typedef LaDa::Vff::VABase<LaDa::Vff::Layered> t_Vff;
# define __PROGNAME__ "Valence Force Field Functional for Epitaxial Zinc-Blende Structures"
#else
# include <vff/functional.h>
# include <vff/va.h>
  typedef LaDa::Vff::VABase<LaDa::Vff::Functional> t_Vff;
# define __PROGNAME__ "Valence Force Field Functional for Zinc-Blende Structures"
#endif

#ifdef LADA_MPI
# include <boost/mpi/environment.hpp>
#endif
    
int main(int argc, char *argv[]) 
{  
  namespace fs = boost::filesystem;
  LADA_MPI_START
  LADA_TRY_BEGIN

  __BPO_START__;
  __BPO_HIDDEN__;
# ifdef _LAYERED
    __BPO_SPECIFICS__("Layered Vff Options")
      ("direction", po::value<std::string>(), 
       "If growth direction is NOT the first cell vector/column, then set it here." );
# endif
  po::options_description all; 
  all.add(generic);
# ifdef _LAYERED
    all.add( specific );
# endif
  po::options_description allnhidden;
  allnhidden.add(all).add(hidden);
  po::positional_options_description p; 
  p.add("input", 1); 
  __BPO_MAP__
  __BPO_HELP_N_VERSION__


  fs::path input( vm["input"].as< std::string >() );
  LADA_DO_NASSERT( not ( fs::is_regular( input ) or fs::is_symlink( input ) ),
              input << " is a not a valid file.\n" );
# ifdef _LAYERED
    LaDa::math::rVector3d direction;
    const bool setdirection = ( vm.count("direction") != 0 );
    if( setdirection )
    {
      boost::tuple<LaDa::types::t_real, LaDa::types::t_real, LaDa::types::t_real> vec;
      LaDa::opt::tuples::read<LaDa::types::t_real>( vm["direction"].as<std::string>(), vec );
      direction[0] = vec.get<0>();
      direction[1] = vec.get<1>();
      direction[2] = vec.get<2>();
    }
# endif
    
  boost::shared_ptr<LaDa::Crystal::Lattice>  
    lattice( LaDa::Crystal::read_lattice( input ) );
  LaDa::Crystal::Structure::lattice = lattice.get();

  TiXmlElement *child;
  LaDa::math::rVector3d vec;
  LaDa::Crystal::Structure structure;
    
  TiXmlDocument doc( input.string().c_str() );
  LADA_DO_NASSERT( not doc.LoadFile(),
              "Error while opening input file.\n";)
  TiXmlHandle handle( &doc );
  
  // loads structure
  child = handle.FirstChild( "Job" ).FirstChild( "Structure" ).Element();
  LADA_DO_NASSERT( not child, "No <Structure> tags were found in input file.\n")
  for (; child ; child = child->NextSiblingElement("Structure") )
  {
    LADA_DO_NASSERT( not structure.Load(*child), 
                "Error while reading Structure description from input.\n")
  
    // loads vff functional
    TiXmlElement *vff_xml = handle.FirstChild( "Job" ).Element();
    LADA_DO_NASSERT( not vff_xml, "Could not find <Job> tag in input.\n")
    t_Vff vff(structure);
    LADA_MPI_CODE( vff.set_mpi( world.get() ); )
    LADA_DO_NASSERT( not vff.Load(*vff_xml),
                "Could not Load Valence Force Field Functional from input.\n")
    if( not vff.init(true) ) continue;
#   ifdef _LAYERED
      if( setdirection ) vff.Vff().set_direction( direction );
#   endif
    
    structure.energy = vff.evaluate() / 16.0217733;

#   ifdef LADA_MPI
      if( world->rank() != 0 ) continue;
#   endif

    const LaDa::math::rMatrix3d stress = vff.Vff().get_stress();
    std::cout << std::fixed << std::setprecision(12) 
              << "Energy [eV]: " << std::setw(18) << structure.energy << std::endl
              << std::setprecision(5)
              << "Stress Tensor: " << std::endl 
              << std::setw(12) << stress(0,0) << " " << std::setw(12)
                               << stress(1,0) << " " << std::setw(12)
                               << stress(2,0) << std::endl
              << std::setw(12) << stress(0,1) << " " << std::setw(12)
                               << stress(1,1) << " " << std::setw(12)
                               << stress(2,1) << std::endl
              << std::setw(12) << stress(0,2) << " " << std::setw(12)
                               << stress(1,2) << " " << std::setw(12)
                               << stress(2,2) 
              << std::endl << std::endl << std::endl;
    vff.print_escan_input( "atomic.config" );
  
  
    LaDa::Crystal::Fourier( structure.atoms.begin(),  structure.atoms.end(),
                       structure.k_vecs.begin(), structure.k_vecs.end() );
    structure.print_out(std::cout );
  }

  return 0;
  __BPO_CATCH__( LADA_MPI_CODE( MPI_Finalize() ) )
}

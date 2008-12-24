//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/program_options.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>
#include <stdexcept>

#include <tinyxml/tinyxml.h>

#include <crystal/lattice.h>
#include <crystal/structure.h>
#include <revision.h>
#include <opt/debug.h>
#include <opt/bpo_macros.h>

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

#include <print/manip.h>

#ifdef _MPI
# include <boost/mpi/environment.hpp>
#endif
    
int main(int argc, char *argv[]) 
{  
  namespace fs = boost::filesystem;
  __MPI_START__
  __TRYBEGIN

  __BPO_START__;
  __BPO_HIDDEN__;
  po::options_description all; 
  all.add(generic);
  po::options_description allnhidden;
  allnhidden.add(all).add(hidden);
  po::positional_options_description p; 
  p.add("input", 1); 
  __BPO_MAP__
  __BPO_HELP_N_VERSION__


  fs::path input( vm["input"].as< std::string >() );
  __DOASSERT( not ( fs::is_regular( input ) or fs::is_symlink( input ) ),
              input << " is a not a valid file.\n" );
    
  boost::shared_ptr<LaDa::Crystal::Lattice>  
    lattice( LaDa::Crystal::read_lattice( input ) );
  LaDa::Crystal::Structure::lattice = lattice.get();

  TiXmlElement *child;
  LaDa::atat::rVector3d vec;
  LaDa::Crystal::Structure structure;
    
  TiXmlDocument doc( input.string().c_str() );
  __DOASSERT( not doc.LoadFile(),
              "Error while opening input file.\n";)
  TiXmlHandle handle( &doc );
  
  // loads structure
  child = handle.FirstChild( "Job" ).FirstChild( "Structure" ).Element();
  __DOASSERT( not child, "No <Structure> tags were found in input file.\n")
  for (; child ; child = child->NextSiblingElement("Structure") )
  {
    __DOASSERT( not structure.Load(*child), 
                "Error while reading Structure description from input.\n")
  
    // loads vff functional
    TiXmlElement *vff_xml = handle.FirstChild( "Job" ).Element();
    __DOASSERT( not vff_xml, "Could not find <Job> tag in input.\n")
    t_Vff vff(structure);
    __DOASSERT( not vff.Load(*vff_xml),
                "Could not Load Valence Force Field Functional from input.\n")
    if( not vff.init(true) ) continue;
    
    structure.energy = vff.evaluate() / 16.0217733;
    const LaDa::atat::rMatrix3d stress = vff.Vff().get_stress();
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
  __BPO_CATCH__( __MPICODE( MPI_Finalize() ) )
}

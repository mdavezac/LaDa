//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/scoped_ptr.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#ifdef _MPI
# include <boost/mpi/environment.hpp>
#endif

#include <tinyxml/tinyxml.h>

#include <print/manip.h>
#include <revision.h>
#include <opt/bpo_macros.h>
#include <opt/initial_path.h>
#include <opt/tinyxml.h>

#include "va.h"
// #include "dipole_elements.h"
#ifdef _LAYERED_
# include <vff/layered.h>
  typedef LaDa::Pescan::VirtualAtom< LaDa::Vff::Layered > t_Pescan;
  typedef LaDa::Vff::VABase<LaDa::Vff::Layered> t_Vff;
# define __PROGNAME__ "Empirical Pseudo-Potential Functional for layered Zinc-Blend Structures"

#else
  typedef LaDa::Pescan::VirtualAtom< LaDa::Vff::Functional > t_Pescan;
  typedef LaDa::Vff::VABase<LaDa::Vff::Functional> t_Vff;
# define __PROGNAME__ "Empirical Pseudo-Potential Functional for Zinc-Blend Structures"
#endif

// performs computation.
struct Eval
{
  LaDa::Crystal::Structure &structure;
  t_Pescan pescan;
  bool doallelectron;
  bool do_evaluate;
// bool compute_dipoles;
// LaDa::types::t_real degeneracy;
  
  Eval( LaDa::Crystal :: Structure & _str ) : structure( _str ), pescan( _str ) {}

  bool operator()( const TiXmlElement &_node )
  {
    LaDa::types::t_real dip;
    if( not structure.Load(_node) ) return false;
    __ROOTCODE
    (
      (*LaDa::mpi::main),
      std::cout << "Successfuly read structure input from file:\n";
    )
    structure.set_site_indices();
    if( not pescan.init( true ) )
    {
      std::cerr << "Error while initializing vff functional. Will skip current structure.\n";
      return false;
    }

    LaDa::Pescan::BandGap& bandgap = pescan.BandGap();
    bandgap.escan.rspace_output = LaDa::Pescan::Interface::Escan::WFN_AFTER_CALL;
    if( doallelectron ) bandgap.set_method( LaDa::Pescan::Interface::ALL_ELECTRON );
    if( do_evaluate ) structure.energy = pescan.evaluate();
//   if( compute_dipoles )
//     dip = LaDa::Pescan::oscillator_strength( bandgap, structure, degeneracy, true );

    LaDa::Crystal::Fourier( structure.atoms.begin(), structure.atoms.end(),
                      structure.k_vecs.begin(), structure.k_vecs.end() );

    __ROOTCODE
    (
      (*LaDa::mpi::main),
      std::cout << structure 
                << "\n\nVBM: " << bandgap.bands.vbm
                << " -- CBM: " << bandgap.bands.cbm
                << "    ---    Band Gap: " << bandgap.bands.gap() << std::endl;
    // if( not compute_dipoles ) return true; 
    // std::cout << "oscillator_strength: " << dip << "\n";
    )
    return true;
  }
};

int main(int argc, char *argv[]) 
{
  LaDa::opt::InitialPath::init();
  namespace fs = boost::filesystem;

  __MPI_START__
  __TRYBEGIN 

  
  __BPO_START__;
  __BPO_HIDDEN__;
  __BPO_SPECIFICS__( "Escan Specific Options" )
#     ifndef _EMASS
//       ( "dipoles,d", "Compute Dipole Moments" )
//       ( "degeneracy", po::value<LaDa::types::t_real>()->default_value(LaDa::types::tolerance),
//         "Allowable band degeneracy when computing dipole moments" )
        ( "diag", "Full diagonalization." )
#     endif
//     ("check,c", po::value<std::string>(), "GA output filename." ) 
      ("donteval", "Whether to perform evaluation." );
  __BPO_GENERATE__()
  __BPO_MAP__ 
  __BPO_HELP_N_VERSION__

  const fs::path filename( vm["input"].as< std::string >() );

  LaDa::Crystal::Structure structure;
  Eval evaluate(structure);
  evaluate.do_evaluate = vm.count( "donteval" ) == 0;
// evaluate.compute_dipoles = vm.count("dipoles") != 0;
// evaluate.degeneracy = vm["degeneracy"].as<LaDa::types::t_real>();
  evaluate.doallelectron = vm.count("diag") != 0;

  __ROOTCODE
  (
    (*LaDa::mpi::main),
    std::cout << "Input filename: " << filename
//             << "\nWill " << ( evaluate.compute_dipoles ? " ": "not " ) 
//             << "compute dipole moments."
              << "\n";
    if( not evaluate.do_evaluate ) 
      std::cout << "Will not perform evaluation.\n";
    std::cout << "\n\n";
  )


  __DOASSERT( not fs::exists( filename ), filename << " does not exist.\n" );
  __DOASSERT( not ( fs::is_regular( filename ) or fs::is_symlink( filename ) ),
              filename << " is a not a valid file.\n" );

  std::string input_filestream;
  LaDa::opt::read_xmlfile( filename, input_filestream );

  boost::shared_ptr< LaDa::Crystal::Lattice >
    lattice( LaDa::Crystal::read_lattice( input_filestream ) );
  LaDa::Crystal::Structure::lattice = lattice.get();


  TiXmlDocument doc;
  TiXmlHandle handle( &doc );
  doc.Parse( input_filestream.c_str() );

  TiXmlElement *child = handle.FirstChild( "Job" ).Element();
  __DOASSERT( not child, "Could not find node \"Job\" in " << filename << ".\n" )
  __DOASSERT( not ( child and evaluate.pescan.Load(*child) ),
              "Error while reading pescan from input.\n" )
  __DIAGA( evaluate.pescan.BandGap().set_mpi( LaDa::mpi::main ); )


  TiXmlHandle docHandle( &doc ); 
  child = docHandle.FirstChild( "Restart" ).
                    FirstChild( "Results" ).
                    FirstChild( "optimum" ).
                    FirstChild( "Indivividual" ).
                    FirstChild( "Structure" ).Element();
  child = handle.FirstChild( "Job" ).FirstChild( "Structure" ).Element();
  for(; child; child = child->NextSiblingElement("Structure" ) )
    evaluate( *child );

  __BPO_CATCH__( __MPICODE( MPI_Finalize() ) )
  catch(...)
  {
    std::cerr << "Finished with an error." << std::endl;
    __MPICODE( MPI_Finalize(); )
    return 1;
  }

  return 0;
}




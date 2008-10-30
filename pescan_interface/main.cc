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


#include "va.h"
#include "dipole_elements.h"
typedef Pescan::VirtualAtom< Vff::Functional > t_Pescan;
typedef Vff::VABase<Vff::Functional> t_Vff;
#define __PROGNAME__ "Empirical Pseudo-Potential Functional"

// performs computation.
struct Eval
{
  Crystal::Structure &structure;
  t_Pescan pescan;
  bool doallelectron;
  bool do_evaluate;
  bool compute_dipoles;
  types::t_real degeneracy;
  
  Eval( Crystal :: Structure & _str ) : structure( _str ), pescan( _str ) {}

  bool operator()( const TiXmlElement &_node )
  {
    types::t_real dip;
    if( not structure.Load(_node) ) return false;
    __ROOTCODE
    (
      (*::mpi::main),
      std::cout << "Successfuly read structure input from file:\n";
    )
    structure.set_site_indices();
    if( not pescan.init( true ) )
    {
      std::cerr << "Error while initializing vff functional. Will skip current structure.\n";
      return false;
    }

    Pescan::BandGap& bandgap = pescan.BandGap();
    bandgap.escan.rspace_output = Pescan::Interface::Escan::WFN_AFTER_CALL;
    if( doallelectron ) bandgap.set_method( Pescan::Interface::ALL_ELECTRON );
    if( do_evaluate ) structure.energy = pescan.evaluate();
    if( compute_dipoles ) dip = Pescan::oscillator_strength( bandgap, structure, degeneracy, true );

    Crystal::Fourier( structure.atoms.begin(), structure.atoms.end(),
                      structure.k_vecs.begin(), structure.k_vecs.end() );

    __ROOTCODE
    (
      (*::mpi::main),
      std::cout << structure 
                << "\n\nVBM: " << bandgap.bands.vbm
                << " -- CBM: " << bandgap.bands.cbm
                << "    ---    Band Gap: " << bandgap.bands.gap() << std::endl;
      if( not compute_dipoles ) return true; 
      std::cout << "oscillator_strength: " << dip << "\n";
    )
    return true;
  }
};

int main(int argc, char *argv[]) 
{
  opt::InitialPath::init();
  namespace fs = boost::filesystem;

  __MPI_START__
  __TRYBEGIN 

  
  __BPO_START__;
  __BPO_SPECIFICS__( "Escan Specific Options" )
#     ifndef _EMASS
        ( "dipoles,d", "Compute Dipole Moments" )
        ( "degeneracy", po::value<types::t_real>()->default_value(types::tolerance),
          "Allowable band degeneracy when computing dipole moments" )
        ( "diag", "Full diagonalization." )
#     endif
      ("check,c", po::value<std::string>(), "GA output filename." ) 
      ("donteval", "Whether to perform evaluation." );
  __BPO_GENERATE__()
  __BPO_MAP__ 
  __BPO_HELP_N_VERSION__

  const fs::path filename( vm["input"].as< std::string >() );
  __DOASSERT( not fs::exists( filename ), filename << " does not exist.\n" );
  __DOASSERT( not ( fs::is_regular( filename ) or fs::is_symlink( filename ) ),
              filename << " is a not a valid file.\n" );
  TiXmlDocument doc( filename.string() );
  TiXmlHandle handle( &doc );
  __ASSERT( not doc.LoadFile(), 
               "error while opening input file " 
            << filename << "\n" << doc.ErrorDesc() << "\n" )

  Crystal::Structure structure;
  Eval evaluate(structure);
  bool do_check( vm.count("check") > 0 );
  fs::path checkfilename;
  if( do_check ) checkfilename = vm["check"].as< std::string >();
  else
  { 
    TiXmlElement *child = handle.FirstChild( "Job" )
                                .FirstChild( "GA" )
                                .FirstChild( "Filenames" )
                                .Element();
    for( ; child; child = child->NextSiblingElement("Filenames") )
      if ( child->Attribute("save")  )
      {
        checkfilename = child->Attribute("save");
        do_check = true;
        __DOASSERT( not ( fs::is_regular( checkfilename ) or fs::is_symlink( checkfilename ) ),
                    checkfilename << " is a not a valid file.\n" );
        break;
      }
  }
  evaluate.do_evaluate = vm.count( "donteval" ) == 0;
  evaluate.compute_dipoles = vm.count("dipoles") != 0;
  evaluate.degeneracy = vm["degeneracy"].as<types::t_real>();
  evaluate.doallelectron = vm.count("diag") != 0;

  __ROOTCODE
  (
    (*::mpi::main),
    std::cout << "Input filename: " << filename
              << "\nWill " << ( evaluate.compute_dipoles ? " ": "not " ) 
              << "compute dipole moments."
              << "\n";
    if( do_check ) std::cout << "Will check GA output structures from file "
                             << checkfilename << "\n";
    if( not evaluate.do_evaluate ) 
      std::cout << "Will not perform evaluation.\n";
    std::cout << "\n\n";
  )



  boost::shared_ptr< Crystal::Lattice >
    lattice( Crystal::read_lattice( filename, "" ) );
  Crystal::Structure::lattice = lattice.get();

  TiXmlElement *child = handle.FirstChild( "Job" ).Element();
  __DOASSERT( not evaluate.pescan.Load(*child),
              "Error while reading pescan from input.\n" )
  __DIAGA( evaluate.pescan.BandGap().set_mpi( mpi::main ); )

  if ( do_check and checkfilename != filename )
    __DOASSERT( not doc.LoadFile(checkfilename.string()),
                "Could not loadfile " << checkfilename << "\n" )


  bool do_check_results( false );
  TiXmlHandle docHandle( &doc ); 
  child = docHandle.FirstChild( "Restart" ).
                    FirstChild( "Results" ).
                    FirstChild( "optimum" ).
                    FirstChild( "Indivividual" ).
                    FirstChild( "Structure" ).Element();
  if( child and evaluate( *child ) ) do_check_results = true;

  child = docHandle.FirstChild( "Restart" ).
                    FirstChild( "Results" ).
                    FirstChild( "Individual").Element();
  for(; child; child = child->NextSiblingElement("Individual" ) )
  {
    if ( not child->FirstChildElement("Structure") ) continue;
    if ( not evaluate( *child ) ) do_check_results = true;
  }
  if( do_check_results ) return 0;
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




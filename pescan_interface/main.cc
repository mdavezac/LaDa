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

#ifdef _EMASS

  #include <vff/layered.h>
  #include <vff/va.h>
  #include "emass.h"
  typedef Pescan::eMassSL t_Pescan;
  typedef Vff::VABase<Vff::Functional> t_Vff;
# define __PROGNAME__ "Empirical Pseudo-Potential Functional for 100 Effective Masses in GaInSb"
  inline void operator<<( t_Pescan &_pescan, const t_Vff &_vff )
  {
    // creates an mpi aware file name for atomic configurations
    std::ostringstream  sstr;
    sstr << "atom_config";
#   ifdef _MPI
      sstr << "." << mpi::main.rank();
#   endif
    // prints atomic configurations
    _vff.print_escan_input(sstr.str());
    // tells bandgap where to find atomic configurations
    _pescan.set_atom_input(sstr.str());
  }

#else

  #include "va.h"
  #include "dipole_elements.h"
  typedef Pescan::VirtualAtom< Vff::Functional > t_Pescan;
  typedef Vff::VABase<Vff::Functional> t_Vff;
# define __PROGNAME__ "Empirical Pseudo-Potential Functional"

#endif



bool evaluate( const TiXmlElement &_node,
               Crystal::Structure &_structure,
               t_Pescan &_pescan, t_Vff &_vff,
               bool _doeval )
{
  if( not _structure.Load(_node) ) return false;
  __ROOTCODE
  (
    (*::mpi::main),
    std::cout << "Successfuly read structure input from file " << std::endl;
  )
  _structure.set_site_indices();

  if( not _vff.init(true) )
  {
    std::cerr << "Error while initializing vff\n" 
              << "Will skip current structure" << std::endl;
    return false;
  }

#ifndef _EMASS
   Pescan::BandGap& bandgap = _pescan.BandGap();
   bandgap.escan.rspace_output = Pescan::Interface::Escan::WFN_AFTER_CALL;
   bandgap.set_method( Pescan::Interface::ALL_ELECTRON );
   if( _doeval ) _structure.energy = _pescan.evaluate();
#else
   _vff.evaluate();
   _pescan << _vff;
   if( _doeval ) _structure.energy = _pescan( _structure );
#endif

  Crystal::Fourier( _structure.atoms.begin(), _structure.atoms.end(),
                     _structure.k_vecs.begin(), _structure.k_vecs.end() );

  return true;
}


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
  bool do_evaluate( vm.count( "donteval" ) == 0 );
# if !defined( _EMASS )
    const bool compute_dipoles( vm.count("dipoles") != 0 );
    const types::t_real degeneracy( vm["degeneracy"].as<types::t_real>() );
# endif

  __ROOTCODE
  (
    (*::mpi::main),
    std::cout << "Input filename: " << filename
#   ifndef _EMASS_ 
              << "\nWill " << ( compute_dipoles ? " ": "not " ) 
              << "compute dipole moments."
#   endif
              << "\n";
    if( do_check ) std::cout << "Will check GA output structures from file "
                             << checkfilename << "\n";
    if( not do_evaluate ) 
      std::cout << "Will not perform evaluation.\n";
    std::cout << "\n\n";
  )



  boost::shared_ptr< Crystal::Lattice >
    lattice( Crystal::read_lattice( filename, "" ) );
  Crystal::Structure::lattice = lattice.get();

  Crystal::Structure structure;
# ifdef _EMASS
    t_Pescan pescan(structure);
    t_Vff vff( structure );
    child = handle.FirstChild( "Job" ).Element();
    if ( not vff.Load(*child) )
    {
      std::cerr << "Error while reading vff from input" << std::endl;
      return false;
    }
# else
    t_Pescan pescan(structure);
    t_Vff &vff = pescan.Vff();
    __DIAGA( pescan.BandGap().set_mpi( mpi::main ); )
# endif
  TiXmlElement *child = handle.FirstChild( "Job" ).Element();
  __DOASSERT( not pescan.Load(*child),
              "Error while reading pescan from input.\n" )

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
  if( child and evaluate( *child, structure, pescan, vff, do_evaluate ) )
  {
    do_check_results = true;
    __ROOTCODE
    (
      (*::mpi::main),
      std::cout << "\n\n\nChecking Optimum\n";
      structure.print_out( std::cout ); 
      std::cout << "\n\n";
      structure.print_xcrysden( std::cout );
#     ifndef _EMASS
        const Pescan::BandGap& bandgap = (const Pescan::BandGap&) pescan;
        std::cout << "\nVBM: " << bandgap.bands.vbm
                  << " -- CBM:" << bandgap.bands.cbm
                  << "    ---    Band Gap: " << bandgap.bands.gap() << std::endl;
      if( compute_dipoles )
        std::cout << "oscillator_strength: "
                  << Pescan::oscillator_strength( bandgap, structure, degeneracy, true )
                  << "\n";
      else std::cout << "no dipole computation\n";
#     else
        std::cout << "Emass tensor:\n" << pescan.tensor;
#     endif
    )
  }

  child = docHandle.FirstChild( "Restart" ).
                    FirstChild( "Results" ).
                    FirstChild( "Individual").Element();
  for(; child; child = child->NextSiblingElement("Individual" ) )
  {
    if ( not child->FirstChildElement("Structure") ) continue;
    if( not evaluate( *child, structure, pescan, vff, do_evaluate ) ) continue;
    do_check_results = true;
    
    __ROOTCODE
    (
      (*::mpi::main),
      std::cout << "\n\n\nNew Structure\n";
      structure.print_out( std::cout ); 
      std::cout << "\n\n";
      structure.print_xcrysden( std::cout );
      if ( not do_evaluate ) continue;
#     ifndef _EMASS
        const Pescan::BandGap& bandgap = (const Pescan::BandGap&) pescan;
        std::cout << "\nVBM: " << bandgap.bands.vbm
                  << " -- CBM:" << bandgap.bands.cbm
                  << "    ---    Band Gap: " << bandgap.bands.gap() << std::endl;
        if( compute_dipoles )
          std::cout << "oscillator_strength: "
                    << Pescan::oscillator_strength( bandgap, structure, degeneracy, true )
                    << "\n";
        else std::cout << "no dipole computation\n";
#     else
        std::cout << "Emass tensor:\n" << pescan.tensor;
#     endif
    )
  }
  if ( do_check_results ) return 0;
  child = handle.FirstChild( "Job" ).FirstChild( "Structure" ).Element();
  for(; child; child = child->NextSiblingElement("Structure" ) )
  {
    if( not evaluate( *child, structure, pescan, vff, do_evaluate ) ) continue;
    
//   __ROOTCODE
//   (
      (*::mpi::main),
      std::cout << "\n\n\nNew Structure\n";
      structure.print_out( std::cout ); 
      std::cout << "\n\n";
      structure.print_xcrysden( std::cout );
//     if ( not do_evaluate ) continue;
#     ifndef _EMASS
//       const Pescan::BandGap& bandgap = (const Pescan::BandGap&) pescan;
        Pescan::BandGap& bandgap = (Pescan::BandGap&) pescan;
        bandgap.eigenvalues.clear();
        bandgap.eigenvalues.push_back( -0.70901723673688466e0 );
        bandgap.eigenvalues.push_back( -0.68309604472960483e0 );
        bandgap.eigenvalues.push_back( -0.67281325704593919e0 );
        bandgap.eigenvalues.push_back( -0.62756363017878858e0 );
        bandgap.eigenvalues.push_back( -0.60229511248203260e0 );
        bandgap.eigenvalues.push_back( -0.58847846945940119e0 );
        bandgap.eigenvalues.push_back( -0.58701852305360425e0 );
        bandgap.eigenvalues.push_back( -0.57925758381173897e0 );
        bandgap.eigenvalues.push_back( -0.44305658641816198e0 );
        bandgap.eigenvalues.push_back( -0.43648201353915156e0 );
        bandgap.eigenvalues.push_back( -0.43393450737341621e0 );
        bandgap.eigenvalues.push_back( -0.43049651938302269e0 );
        bandgap.eigenvalues.push_back( -0.42973720813485394e0 );
        bandgap.eigenvalues.push_back( -0.42935929108758714e0 );
        bandgap.eigenvalues.push_back( -0.42271399633501422e0 );
        bandgap.eigenvalues.push_back( -0.30217478051297092e0 );
        bandgap.eigenvalues.push_back( -0.29759083131265834e0 );
        bandgap.eigenvalues.push_back( -0.29208104213309954e0 );
        bandgap.eigenvalues.push_back( -0.28901572497824213e0 );
        bandgap.eigenvalues.push_back( -0.28484493288616974e0 );
        bandgap.eigenvalues.push_back( -0.28162144758675789e0 );
        bandgap.eigenvalues.push_back( -0.25065308270200265e0 );
        bandgap.eigenvalues.push_back( -0.24758386957743173e0 );
        bandgap.eigenvalues.push_back( -0.24267801736958261e0 );
        bandgap.eigenvalues.push_back( -0.23886634318751651e0 );
        bandgap.eigenvalues.push_back( -0.23290426403269765e0 );
        bandgap.eigenvalues.push_back( -0.22733513805707850e0 );
        bandgap.eigenvalues.push_back( -0.22039331158218275e0 );
        bandgap.eigenvalues.push_back( -0.21810395591884335e0 );
        bandgap.eigenvalues.push_back( -0.20972454395681303e0 );
        bandgap.eigenvalues.push_back( -0.18884797109807500e0 );
        bandgap.eigenvalues.push_back( -0.18294040863380412e0 );
        bandgap.eigenvalues.push_back( -0.18103982066164309e0 );
        bandgap.bands.vbm = -0.18294040863380412e0;
        bandgap.bands.cbm = -0.18103982066164309e0;
        std::cout << "\nVBM: " << bandgap.bands.vbm
                  << " -- CBM:" << bandgap.bands.cbm
                  << "    ---    Band Gap: " << bandgap.bands.gap() << std::endl;

        std::cout << compute_dipoles << "\n";
        if( compute_dipoles )
          std::cout << "oscillator_strength: "
                    << Pescan::oscillator_strength( bandgap, structure, degeneracy, true )
                    << "\n";
#     else
        std::cout << "Emass tensor:\n" << pescan.tensor;
#     endif
  // )
  }

  __BPO_CATCH__

  return 0;
}

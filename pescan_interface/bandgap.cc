//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sstream>
#include <algorithm>

#include <physics/physics.h>
#include <print/stdout.h>
#include <opt/fuzzy.h>
#include <print/manip.h>
#include <opt/debug.h>
#include <opt/initial_path.h>

#include "bandgap.h"

namespace Pescan
{

#ifdef _NOLAUNCH
  //! Fake and fast functional, when compiled with --enable-nolaunch
  void nolaunch_functional( const Crystal::Structure &_str, Bands &bands );
#endif

  types::t_real BandGap::folded_spectrum(const Crystal::Structure &_str)
  {
    const t_Path olddirname( dirname );
    create_directory();

    // First compute CBM
    computation = CBM;
    escan.Eref = Eref.cbm;
    if ( not do_destroy_dir ) 
    {
      dirname = olddirname / "cbm";
      create_directory();
    }
    create_potential();
    __TRYCODE( launch_pescan();
               read_and_throw();,
               "Error while computing CBM.\n") 
    bands.cbm = find_closest_eig( Eref.cbm );

    // Then computes VBM
    computation = VBM;
    escan.Eref = Eref.vbm;
    if ( not do_destroy_dir ) 
    {
      dirname = olddirname / "vbm";
      create_directory();
      create_potential();
    }
    __TRYCODE( launch_pescan();
               read_and_throw();,
               "Error while computing VBM.\n") 
    bands.vbm = find_closest_eig( Eref.cbm );

#ifndef _NOLAUNCH
    if( do_correct and bands.gap() < metallicity ) correct( olddirname );
    if( bands.gap() < metallicity )
    {
      escan.method = ALL_ELECTRON;
      Print :: out << "Failed to find non-metallic band gap.\n"
                   << "Will try an all-electron calculation." << Print::endl;
      all_electron( _str );
    }
#else
    nolaunch_functional( _str, bands );
#endif // _NOLAUNCH

    destroy_directory_();
    dirname = olddirname;

    return bands.gap();
  }

  void BandGap :: correct(const t_Path & _dir) 
  {
    Escan saved_state = escan;
    Bands keeprefs = Eref;
    escan.wavefunction_in = escan.wavefunction_out;
    escan.read_in.resize( escan.nbstates );
    std::vector<types::t_unsigned> :: iterator i_r = escan.read_in.begin();
    std::vector<types::t_unsigned> :: iterator i_r_end = escan.read_in.end();
    for(types::t_unsigned u=1; i_r != i_r_end; ++i_r, ++u ) *i_r = u;

    types::t_unsigned i = 0;
    do 
    {
      Print::out __DODEBUGCODE(<< __SPOT_ERROR)
                 << " Found metallic band gap: " << bands.gap() << "\n" 
                 << "    VBM=" << bands.vbm << "  CBM=" << bands.cbm << "\n" 
                 << " Will Try and modify references.\n";

      if( std::abs( bands.vbm - Eref.vbm ) > std::abs( bands.cbm - Eref.cbm ) )
      {
        Eref.vbm -= inc_dec;
        computation = VBM;
        escan.Eref = Eref.vbm;
        Print::out << "Modifying VBM reference: " << Eref.vbm
                   << Print::endl << Print::flush;
      }
      else
      {
        Eref.cbm += inc_dec;
        computation = CBM;
        escan.Eref = Eref.cbm;
        Print::out << "Modifying CBM reference: " << Eref.cbm
                   << Print::endl << Print::flush;
      }
      if ( not do_destroy_dir )
        dirname = _dir / ( computation == CBM ?  "/cbm": "/vbm" );
      do_correct = false;
      
      __TRYCODE( launch_pescan();
                 read_and_throw();,
                    "Error while computing"
                 << ( computation == CBM ?  "/cbm": "/vbm" )
                 << "\n" )
      do_correct = true;
      computation == VBM ?  bands.vbm = find_closest_eig( Eref.vbm ):
                            bands.cbm = find_closest_eig( Eref.cbm );
      Print::out << "After Correction: VBM=" << bands.vbm
                 << " CBM=" << bands.cbm << Print::endl;
      ++i;
      if ( not do_destroy_dir ) dirname = _dir;
    } while ( bands.gap() < metallicity and i < 5 );

    escan = saved_state;
    Eref  = keeprefs;
  }
  
  types::t_real BandGap::all_electron( const Crystal::Structure &_str ) 
  {
    Crystal::Structure::t_Atoms::const_iterator i_atom = _str.atoms.begin();
    Crystal::Structure::t_Atoms::const_iterator i_atom_end = _str.atoms.end();
    escan.nbstates = 0; 
    for(; i_atom != i_atom_end; ++i_atom)
    {
      Crystal::StrAtom atom; 
      _str.lattice->convert_Atom_to_StrAtom( *i_atom, atom );
      escan.nbstates += Physics::Atomic::Charge( atom.type );
    }
    escan.nbstates += 2;
    if (    escan.potential != Escan::SPINORBIT
         or atat::norm2(escan.kpoint) < types::tolerance )
      escan.nbstates /= 2;

    __TRYCODE( Interface::operator()();,
               "All-electron calculation failed.\n" )

#ifndef _NOLAUNCH
    bands.cbm = eigenvalues[ escan.nbstates - 1 ];
    bands.vbm = eigenvalues[ escan.nbstates - 2 ];
    std::cout << "BandGap: " << bands.gap() << " = " 
              << bands.cbm << " - " << bands.vbm << std::endl;
#else
    nolaunch_functional( _str, bands );
#endif // _NOLAUNCH

    destroy_directory_();

    return bands.gap();
  }

  bool BandGap :: Load( const TiXmlElement &_node )
  {
    const TiXmlElement *parent = find_node( _node );
    __DOASSERT( not parent, 
                "Could not find <Functional type=\"escan\"> in input\n"; )
    __TRYCODE( __DOASSERT( not Interface :: Load_( *parent ), "" ),
               "Encountered error while loading bandgap interface:\n" )
    
    t_method m = escan.method;
    escan.method = ALL_ELECTRON;
    const TiXmlElement *child = parent->FirstChildElement("Reference");
    if( not child ) child = parent->FirstChildElement("References"); 
    if( not child ) return true;
    if(     (not child->Attribute("VBM") )
        and (not child->Attribute("CBM") ) ) return true;
    if(     (not child->Attribute("VBM") )
        or (not child->Attribute("CBM") ) )
    {
      std::cerr __DODEBUGCODE(<< __SPOT_ERROR)
                << "Found only one VBM/CBM reference when expecting two in input \n"
                << "Will continue with all electron calculation." << std::endl;
      return true;
    }

    double d;
    child->Attribute( "CBM", &d );
    Eref.cbm = (types::t_real) d;
    child->Attribute( "VBM", &d );
    Eref.vbm = (types::t_real) d;

    if( Eref.gap() < 0.0 )  
    {
      std::cerr __DODEBUGCODE(<< __SPOT_ERROR)
                << "Found strange VBM/CBM references in input \n"
                << "VBM = " <<  Eref.vbm << " > CBM = " << Eref.cbm << "\n" 
                << "Will continue with all electron calculation." << std::endl;
      return true;
    }

    escan.method = m;
    return true;
  }

  void BandGap :: read_and_throw()
  {
    try { read_result(); }
    catch( std::exception &_e )
    { 
      const t_Path filepath( opt::InitialPath::path() / dirname / escan.output.filename() );
      std::ifstream file;
      file.open( filepath.string().c_str(), std::ios_base::in ); 
      
      __DOASSERT( not file.is_open(),
                    "Could not open escan output file "
                 << filepath << "\n" );
      __DOASSERT( file.bad(),
                    "Error while opening escan output file "
                 << filepath << "\n" );
      char cline[256];
      std::string line("");
      std::cout << "\n\n--- ESCAN OUTPUT ---\n";
      while( not file.eof() )
      {
        file.getline( cline, 256 );
        line = cline;
        std::cout << line << "\n";
      }
      std::cout  << "\n\n--- ESCAN OUTPUT ---\n" << std::endl;
      __THROW_ERROR( _e.what() )
    }
  }

#ifdef _NOLAUNCH
  void nolaunch_functional( const Crystal::Structure &_str, Bands &bands )
  {
    Crystal::Structure::t_kAtoms::const_iterator i_k = _str.k_vecs.begin();
    Crystal::Structure::t_kAtoms::const_iterator i_k_end = _str.k_vecs.end();
    bands.vbm = 0.0; bands.cbm = 0.0;
    bool which = true, sign = true;
    for(; i_k != i_k_end; ++i_k, which = not which )
    {
      types::t_real a0 = std::abs( i_k->type ) * i_k->pos(0);
      types::t_real a1 = std::abs( i_k->type ) * i_k->pos(1);
      types::t_real a2 = std::abs( i_k->type ) * i_k->pos(2);

      bands.vbm +=  a0 * a0 / 5.0 - a1 / 15.0 - a2 * a1;
      bands.cbm += 7.0 * cos( a0 * a1 / 7.0 + a2  ); 
    }

    bands.cbm += bands.vbm;

    if ( Fuzzy::gt( bands.cbm, bands.vbm ) ) return;
    std::swap( bands.cbm, bands.vbm );
    if ( Fuzzy::eq(bands.vbm, bands.cbm ) ) bands.cbm += 0.1;
  }
#endif

}

//
//  Version: $Id$
//
#include <algorithm>

#include <physics/physics.h>
#include <print/stdout.h>
#include <opt/fuzzy.h>

#include "bandgap.h"

namespace Pescan
{
#ifdef _NOLAUNCH
  //! Fake and fast functional, when compiled with --enable-nolaunch
  void nolaunch_functional( const Ising_CE::Structure &_str, Bands &bands );
#endif

#ifdef _NOLAUNCH
  types::t_real BandGap::folded_spectrum(const Ising_CE::Structure &_str)
#else
  types::t_real BandGap::folded_spectrum()
#endif
  {
    std::string olddirname = dirname;
    create_directory();

    // First compute CBM
    computation = CBM;
    escan.Eref = Eref.cbm;
    if ( not do_destroy_dir ) 
    {
      dirname = olddirname + ( computation == CBM ?  "/cbm": "/vbm" );
      create_directory();
    }
    create_potential();
    launch_pescan();
    read_result();
    bands.cbm = find_closest_eig( Eref.cbm );
    destroy_directory_();

    // Then computes VBM
    computation = VBM;
    escan.Eref = Eref.vbm;
    if ( not do_destroy_dir ) 
      dirname = olddirname + ( computation == CBM ?  "/cbm": "/vbm" );
    create_directory();
    launch_pescan();
    read_result();
    bands.vbm = find_closest_eig( Eref.cbm );
    destroy_directory_();

#ifndef _NOLAUNCH
    if( do_correct and bands.gap() < metallicity ) correct( olddirname );
#else
    nolaunch_functional( _str, bands );
#endif // _NOLAUNCH

    destroy_directory_();
    dirname = olddirname;

    return bands.gap();
  }

  void BandGap :: correct(const std::string & _dir) 
  {
    Escan saved_state = escan;
    Bands keeprefs = Eref;
    std::string wfn = escan.wavefunction_in;
    escan.wavefunction_in = escan.wavefunction_out;
    escan.read_in.resize( escan.nbstates );
    std::vector<types::t_unsigned> :: iterator i_r = escan.read_in.begin();
    std::vector<types::t_unsigned> :: iterator i_r_end = escan.read_in.end();
    for(types::t_unsigned u=0; i_r != i_r_end; ++i_r, ++u ) *i_r = u;

    do 
    {
#ifdef _DEBUG
      Print::out << __FILE__ << ", line: " << __LINE__ << "\n";
#endif
      Print::out << " Found metallic band gap!! " << bands.gap() << "\n" 
                 << " Will Try and modify references. \n";

      if( std::abs( bands.vbm - Eref.vbm ) > std::abs( bands.cbm - Eref.cbm ) )
      {
        Eref.vbm -= inc_dec;
        computation = VBM;
        escan.Eref = Eref.vbm;
      }
      else
      {
        Eref.cbm += inc_dec;
        computation = CBM;
        escan.Eref = Eref.cbm;
      }
      if ( not do_destroy_dir )
        dirname = _dir + ( computation == CBM ?  "/cbm": "/vbm" );
      launch_pescan();
      read_result();
      computation == VBM ?  bands.vbm = find_closest_eig( Eref.vbm ):
                            bands.cbm = find_closest_eig( Eref.cbm );
    } while ( bands.gap() < metallicity );

    escan = saved_state;
    Eref  = keeprefs;
  }
  
  types::t_real BandGap::all_electron( const Ising_CE::Structure &_str ) 
  {
    Ising_CE::Structure::t_Atoms::const_iterator i_atom = _str.atoms.begin();
    Ising_CE::Structure::t_Atoms::const_iterator i_atom_end = _str.atoms.end();
    escan.nbstates = 0; 
    for(; i_atom != i_atom_end; ++i_atom)
    {
      Ising_CE::StrAtom atom; 
      _str.lattice->convert_Atom_to_StrAtom( *i_atom, atom );
      escan.nbstates += Physics::Atomic::Charge( atom.type );
    }
    if (    escan.potential != Escan::SPINORBIT
         or atat::norm2(escan.kpoint) < types::tolerance )
      escan.nbstates /= 2;
    ++escan.nbstates;

    if ( not Interface::operator()() ) return -1;

#ifndef _NOLAUNCH
    bands.cbm = eigenvalues[ escan.nbstates - 1 ];
    bands.vbm = eigenvalues[ escan.nbstates - 2 ];
#else
    nolaunch_functional( _str, bands );
#endif // _NOLAUNCH

    destroy_directory_();

    return bands.gap();
  }

  bool BandGap :: Load( const TiXmlElement &_node )
  {
    const TiXmlElement *parent = find_node( _node );
    if ( not parent )
    {
      std::cerr << "Could not find <Functional type=\"escan\"> in input" << std::endl;
      return false;
    }
    if ( not Interface :: Load_( *parent ) )
    {
      std::cerr << "Could not load Pescan functional" << std::endl;
      return false;
    }

    
    t_method m = escan.method;
    escan.method = ALL_ELECTRON;
    const TiXmlElement *child = parent->FirstChildElement("Reference");
    if( not child ) return true;
    if(     (not child->Attribute("VBM") )
        and (not child->Attribute("CBM") ) ) return true;
    if(     (not child->Attribute("VBM") )
        or (not child->Attribute("CBM") ) )
    {
      std::cerr << "Found only one VBM/CBM reference when expecting two in input \n"
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
      std::cerr << "Found strange VBM/CBM references in input \n"
                << "VBM = " <<  Eref.vbm << " > CBM = " << Eref.cbm << "\n" 
                << "Will continue with all electron calculation." << std::endl;
      return true;
    }

    escan.method = m;
    return true;
  }


#ifdef _NOLAUNCH
  void nolaunch_functional( const Ising_CE::Structure &_str, Bands &bands )
  {
    Ising_CE::Structure::t_kAtoms::const_iterator i_k = _str.k_vecs.begin();
    Ising_CE::Structure::t_kAtoms::const_iterator i_k_end = _str.k_vecs.end();
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

    if ( Fuzzy::ge( bands.cbm, bands.vbm ) ) return;
    std::swap( bands.cbm, bands.vbm );
    if ( Fuzzy::eq(bands.vbm, bands.cbm ) ) bands.cbm += 0.1;
  }
#endif

}

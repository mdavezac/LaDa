//
//  Version: $Id$
//
#include <physics/physics.h>
#include <print/stdout.h>

#include "bandgap.h"

namespace Pescan
{
  types::t_real BandGap::folded_spectrum()
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

    if( do_correct and bands.gap() < 0.001 ) correct( olddirname );

    destroy_directory_();
    dirname = olddirname;

    return bands.gap();
  }

  void BandGap :: correct(const std::string & _dir) 
  {
    Bands keeprefs = Eref;
    bool old_input = do_input_wavefunctions;
    do_input_wavefunctions = true;
    std::string wfn = escan.wavefunction_in;
    escan.wavefunction_in = escan.wavefunction_out;
    do 
    {
      Print::out << " Found metallic band gap!! " << bands.gap() << "\n" 
                 << " Will Try and modify references. \n";

      if( std::abs( bands.vbm - Eref.vbm ) > std::abs( bands.cbm - Eref.cbm ) )
      {
        Eref.vbm -= 0.05;
        computation = VBM;
        escan.Eref = Eref.vbm;
      }
      else
      {
        Eref.cbm += 0.05;
        computation = CBM;
        escan.Eref = Eref.cbm;
      }
      if ( not do_destroy_dir )
        dirname = _dir + ( computation == CBM ?  "/cbm": "/vbm" );
      launch_pescan();
      read_result();
      computation == VBM ?  bands.vbm = find_closest_eig( Eref.vbm ):
                            bands.cbm = find_closest_eig( Eref.cbm );
    } while ( bands.gap() < 0.001 );

    do_input_wavefunctions = old_input;
    escan.wavefunction_in = wfn;
    Eref = keeprefs;
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

    bands.cbm = eigenvalues[ escan.nbstates - 1 ];
    bands.vbm = eigenvalues[ escan.nbstates - 2 ];

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

}

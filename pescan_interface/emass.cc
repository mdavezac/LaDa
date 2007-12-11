//
//  Version: $Id$
//
#include <physics/physics.h>
#include <print/stdout.h>

#include "emass.h"

extern "C" void FC_FUNC(dsyev, DSYEV)( char*, char*, int*, types::t_real*, 
                                       int*, types::t_real*, types::t_real*, 
                                       int*, int* );
                                                    
namespace Pescan
{
  const atat::rVector3d eMassSL :: Gamma(0,0,0);
  const atat::rVector3d eMassSL :: Lx(1,0,0);
  const atat::rVector3d eMassSL :: Ly(0,1,0);
  const atat::rVector3d eMassSL :: Lz(0,0,1);
  const atat::rVector3d eMassSL :: Hxy(eMassSL::sqrt2,eMassSL::sqrt2,0);
  const atat::rVector3d eMassSL :: Hxmy(eMassSL::sqrt2,eMassSL::msqrt2,0);
  const atat::rVector3d eMassSL :: Hyz(0,eMassSL::sqrt2,eMassSL::sqrt2);
  const atat::rVector3d eMassSL :: Hymz(0,eMassSL::sqrt2,eMassSL::msqrt2);
  const atat::rVector3d eMassSL :: Hxz(eMassSL::sqrt2,0,eMassSL::sqrt2);
  const atat::rVector3d eMassSL :: Hmxz(eMassSL::msqrt2,0,eMassSL::sqrt2);

  bool eMassSL::operator()( const Ising_CE::Structure &_str )
  {
    escan_calls( _str );

    const types::t_real amp2 =   _str.scale * _str.scale 
                               * Physics::emass("kg") / Physics::hbar("eV*s") 
                               / Physics::hbar("J*s") / 4.0 
                               / Math::pi / Math::pi / amplitude / amplitude;
    inverse(0,0) = ( eig_Lx[0] + eig_Lx[1] - 2.0 * eig_gamma ) * amp2;
    inverse(1,1) = ( eig_Ly[0] + eig_Ly[1] - 2.0 * eig_gamma ) * amp2;
    inverse(2,2) = ( eig_Lz[0] + eig_Lz[1] - 2.0 * eig_gamma ) * amp2;
    inverse(0,1) = inverse(1,0) = (   eig_Hxy[0] - eig_Hxmy[0]
                                    + eig_Hxy[1] - eig_Hxmy[1] ) * amp2;
    inverse(0,2) = inverse(2,0) = (   eig_Hyz[0] - eig_Hymz[0]
                                    + eig_Hyz[1] - eig_Hymz[1] ) * amp2;
    inverse(1,2) = inverse(2,1) = (   eig_Hxz[0] - eig_Hmxz[0]
                                    + eig_Hxz[1] - eig_Hmxz[1] ) * amp2;

    types::t_real work[9];
    types::t_real eigs[3];
    types::t_real inv[9];
    for( types::t_unsigned i=0; i<3; ++i )
      for( types::t_unsigned j=0; j<3; ++j )
      {
        if( opt::Fuzzy<types::t_real>::equal( inverse(i,j), 0.0 ) ) 
          inverse(i,j) = 0.0;
        inv[3*i+j] = inverse(i,j);
      }
    

    int info;
    char job('V'), mattype('U');
    int three(3), nine(9);
    FC_FUNC(dsyev, DSYEV)( &job, &mattype, &three, inv, &three, eigs, work, &nine, &info);
    if( info != 0 ) return false;

    atat::rMatrix3d m, vecs; m.zero();
    m(0,0) = eigs[0];
    m(1,1) = eigs[1];
    m(2,2) = eigs[2];

    for( types::t_unsigned i=0; i<3; ++i )
      for( types::t_unsigned j=0; j<3; ++j )
        vecs(i,j) = inv[3*i+j];

    tensor = vecs * m * (~vecs);
    std::cout << " Inverse mass: \n " << inverse
              << " Inverse mass: \n " << tensor << "\n";

    m.zero();
    m(0,0) = 1.0 / eigs[0];
    m(1,1) = 1.0 / eigs[1];
    m(2,2) = 1.0 / eigs[2];

    tensor = vecs * m * (~vecs);
    std::cout << " mass: \n " << tensor << "\n";
    return true;
  }

  void eMassSL :: escan_calls( const Ising_CE::Structure &_str )
  {
    // Stuff to save
    atat::rVector3d vec = escan.kpoint;
    types::t_unsigned nb = escan.nbstates;
    std::string in = escan.wavefunction_in;
    std::string out = escan.wavefunction_out;

    // Needed only once
    create_directory();
    create_potential();

    // Compute Gamma. This is the one computation which can be all electron
    eig_gamma =   Physics::Hartree("eV")
                * ( escan.method == Interface::FOLDED_SPECTRUM ) ?
                   gamma_folded_spectrum():
                   gamma_all_electron( _str );

    // Then compute other points
    bool m = do_input_wavefunctions;
    do_input_wavefunctions = true;
    escan.nbstates = 2;
    escan.wavefunction_in = escan.wavefunction_out;
    escan.wavefunction_out = "dummy";
   
    other_kpoints( amplitude * Lx, eig_Lx);
    other_kpoints( amplitude * Ly, eig_Ly);
    other_kpoints( amplitude * Lz, eig_Lz);
    other_kpoints( amplitude * Hxy, eig_Hxy);
    other_kpoints( amplitude * Hxmy, eig_Hxmy);
    other_kpoints( amplitude * Hyz, eig_Hyz);
    other_kpoints( amplitude * Hymz, eig_Hymz);
    other_kpoints( amplitude * Hxz, eig_Hxz);
    other_kpoints( amplitude * Hmxz, eig_Hmxz);

    do_input_wavefunctions = m;
    escan.kpoint = vec;
    escan.nbstates = nb;
    escan.wavefunction_in = in;
    escan.wavefunction_out = out;
    destroy_directory_();
  }

  types::t_real eMassSL::gamma_all_electron( const Ising_CE::Structure &_str ) 
  {
    escan.kpoint = Gamma;
    
    Ising_CE::Structure::t_Atoms::const_iterator i_atom = _str.atoms.begin();
    Ising_CE::Structure::t_Atoms::const_iterator i_atom_end = _str.atoms.end();
    escan.nbstates = 2; 
    for(; i_atom != i_atom_end; ++i_atom)
    {
      Ising_CE::StrAtom atom; 
      _str.lattice->convert_Atom_to_StrAtom( *i_atom, atom );
      escan.nbstates += Physics::Atomic::Charge( atom.type );
    }
    if (    escan.potential != Escan::SPINORBIT
         or atat::norm2(escan.kpoint) < types::tolerance )
      escan.nbstates /= 2;

    launch_pescan();
    read_result();

    escan.Eref = eigenvalues.back();
    return escan.Eref;
  }

  bool eMassSL :: Load( const TiXmlElement &_node )
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
    if(     (not child->Attribute("value") )
        and (not child->Attribute("VBM") ) ) return true;

    double d;
    if( child->Attribute( "value" ) ) child->Attribute( "value", &d );
    if( child->Attribute( "CBM" ) ) child->Attribute( "CBM", &d );
    escan.Eref = (types::t_real) d;
    escan.method = m;
    return true;
  }

}

//
//  Version: $Id$
//
// Defines two new classes, 
//   Atom: object containing a position (as atat::rVector3d)
//            and a type (types::t_real)
//   Structure: object containing a cell (as atat::rMatrix3d) 
//                 and an array of Atom (in std::vector container)


#ifndef _ISING_CE_STRUCTURE_H_
#define _ISING_CE_STRUCTURE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include <ostream>
#include <fstream>
#include <string>
#include <complex>
#include <math.h>


#include <tinyxml/tinyxml.h>
#include <exception>

#include "opt/types.h"
#include "atat/vectmac.h"
#include "atat/xtalutil.h"
#include "atat/machdep.h"
#include "atom.h"
#include "lattice.h"

namespace Ising_CE {

  typedef atat::Structure Atat_Structure;
  
  struct Structure
  {
    typedef Atom_Type<types::t_real>  t_Atom;
    typedef std::vector< t_Atom >     t_Atoms;
    typedef CAtom                     t_kAtom;
    typedef std::vector< t_kAtom >    t_kAtoms;

    const static types::t_unsigned FREEZE_NONE;
    const static types::t_unsigned FREEZE_XX;
    const static types::t_unsigned FREEZE_XY;
    const static types::t_unsigned FREEZE_XZ;
    const static types::t_unsigned FREEZE_YY;
    const static types::t_unsigned FREEZE_YZ;
    const static types::t_unsigned FREEZE_ZZ;

    atat::rMatrix3d cell;
    std::vector< Atom_Type<types::t_real> > atoms;
    std::vector< CAtom > k_vecs;
    types::t_int Pi_name;
    types::t_real energy;
    types::t_unsigned freeze;
    types::t_real scale;
    static Lattice *lattice;

    Structure() : Pi_name(0), energy(0), freeze(FREEZE_NONE) {};
    Structure   ( const Atat_Structure &atat ) 
              : Pi_name(0), energy(0), freeze(FREEZE_NONE)
        { convert_from_ATAT( atat ); };
    Structure   ( const TiXmlElement &_element )
              : Pi_name(0), energy(0), freeze(FREEZE_NONE)
      { Load( _element ); };
    Structure   ( const Structure &_str )
              : cell(_str.cell), atoms(_str.atoms), k_vecs(_str.k_vecs),
                Pi_name(_str.Pi_name), energy(_str.energy), freeze(_str.freeze) {}
    ~Structure () {};
    void convert_from_ATAT (const Atat_Structure &atat);
    void operator=( const Atat_Structure &atat )
     { convert_from_ATAT( atat ); };
    void print_out( std::ostream &stream ) const;
    void set_atom_types( const std::vector<types::t_real> &types);
    void get_atom_types( std::vector<types::t_real> &types) const;
    types::t_real get_concentration() const;
    void randomize(types::t_real range, bool centered);
    bool Load( const TiXmlElement &_element );
    void print_xml( TiXmlElement &_node ) const;
    bool operator== (const Structure &_str ) const
    {
      if ( not (cell == _str.cell) )
        return false;
      return atoms == _str.atoms;
    }
    
    template<class t_container >
    void set_kvectors( const t_container &_container )
    {
      typename t_container :: const_iterator i_kvec =  _container.begin();
      typename t_container :: const_iterator i_end =  _container.end();
      CAtom kvec;
      k_vecs.clear();
      k_vecs.reserve( _container.size() );
      for( ; i_kvec != i_end; ++i_kvec ) 
      {
        kvec = (*i_kvec);
        k_vecs.push_back( kvec );
      }
    }

    bool set_site_indices();

    void find_k_vectors();
  };

  void  find_range( const atat::rMatrix3d &A, atat::iVector3d &kvec );
  void refold( atat::rVector3d &vec, const atat::rMatrix3d &lat );
  bool are_equivalent( const atat::rVector3d &_a,
                       const atat::rVector3d &_b,
                       const atat::rMatrix3d &_cell);

  template <class CONTAINER>
  void remove_equivalents( CONTAINER &_cont, const atat::rMatrix3d &_cell)
  {
    typename CONTAINER :: iterator i_vec = _cont.begin();
    typename CONTAINER :: iterator i_end = _cont.end();
    typename CONTAINER :: iterator i_which;

    while( i_vec != i_end )
    {
      i_which = i_vec+1;
      for ( ; i_which != i_end; i_which++ )
        if ( are_equivalent( *i_which, *i_vec, _cell ) )
          break;

      if ( i_which == i_end )
      { 
        ++i_vec;
        continue;
      }
      
      ( atat::norm2( (atat::rVector3d&) *i_vec ) < atat::norm2( (atat::rVector3d&) *i_which ) ) ? 
            _cont.erase(i_which): _cont.erase(i_vec);
      i_vec = _cont.begin();
      i_end = _cont.end();
    }

  }
  template<class T_R_IT, class T_K_IT>
  void fourrier_to_kspace( T_R_IT _rfirst, T_R_IT _rend,
                           T_K_IT _kfirst, T_K_IT _kend ) // sets kvector values from rspace values
  {
    const std::complex<types::t_real> imath(0, -2*3.1415926535897932384626433832795028841971693993751058208);
    
    for (; _kfirst != _kend; ++_kfirst)
    {
      _kfirst->type = std::complex<types::t_real>(0);
      for(T_R_IT i_r( _rfirst ); i_r != _rend; ++i_r )
      {
        _kfirst->type +=   exp( imath * ( i_r->pos[0] * _kfirst->pos[0] +
                                          i_r->pos[1] * _kfirst->pos[1] +
                                          i_r->pos[2] * _kfirst->pos[2] ) )
                         * i_r->type;
      }
    }
  }
  template<class T_R_IT, class T_K_IT, class T_O_IT >
  void fourrier_to_rspace( T_R_IT _rfirst, T_R_IT _rend,
                           T_K_IT _kfirst, T_K_IT _kend,
                           T_O_IT _rout ) // sets rvector values from kspace values
  {
    const std::complex<types::t_real> imath(0, 2*3.1415926535897932384626433832795028841971693993751058208);
    for (; _rfirst != _rend; ++_rfirst, ++_rout)
    {
      *_rout = 0.0;
      for(T_K_IT i_k=_kfirst; i_k != _kend; ++i_k)
      {
        *_rout +=   exp( imath * ( _rfirst->pos[0] * i_k->pos[0] +
                                   _rfirst->pos[1] * i_k->pos[1] +
                                   _rfirst->pos[2] * i_k->pos[2] ) )
                  * i_k->type;
      }
    }
  }

  inline std::ostream& operator<<( std::ostream& _stream, const Ising_CE::Structure& _struc )
    { _struc.print_out(_stream); return _stream; }

} // namespace Ising_CE


#endif

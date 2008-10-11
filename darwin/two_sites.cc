//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <algorithm>

#include <opt/fuzzy.h>
#include <opt/debug.h> 
#include <opt/atat.h>

#include "two_sites.h"

namespace TwoSites
{
  void rearrange_structure( Crystal::Structure &_str)
  {
    if ( not _str.lattice and _str.lattice->sites.size() != 2)
      return;

    std::list< Crystal::Structure::t_Atom > sites0;
    std::list< Crystal::Structure::t_Atom > sites1;
    Crystal::Structure::t_Atoms :: const_iterator i_atom = _str.atoms.begin();
    Crystal::Structure::t_Atoms :: const_iterator i_atom_end = _str.atoms.end();
    for(; i_atom != i_atom_end; ++i_atom )
      ( i_atom->site == 0 ) ?
        sites0.push_back( *i_atom ): sites1.push_back( *i_atom );

    std::list< Crystal::Structure::t_Atom > :: iterator i_0 = sites0.begin();
    std::list< Crystal::Structure::t_Atom > :: iterator i_end = sites0.end();
    std::list< Crystal::Structure::t_Atom > :: iterator i_1;
    atat::rVector3d translation =   _str.lattice->sites[1].pos
                                  - _str.lattice->sites[0].pos; 
    _str.atoms.clear();
    for(; i_0 != i_end; ++i_0 )
    {
      atat::rVector3d atom = i_0->pos + translation;
      i_1 = std::min_element( sites1.begin(), sites1.end(), atat::norm_compare( atom ) );
      _str.atoms.push_back( *i_0 ); 
      _str.atoms.push_back( *i_1 ); 
    }

  }

  //  sets concentration from k-space values.
  //  individual need not be physical (ie S_i=+/-1) when fourier transformed to real-space
  //  a normalization procedure is applied which takes into account:
  //  (i) that this ga run is at set concentration (x =x0, y=y0)
  //  (ii) that x and y are only constrained by load-balancing
  void Concentration :: operator()( Crystal::Structure &_str )
  {
    types::t_complex  *hold = new types::t_complex[ N ];
    __DOASSERT( not hold,  "Memory allocation error\n" )

    // creates individual with unnormalized occupation
    types::t_complex  *i_hold = hold;
    TwoSites::Fourier( _str.atoms.begin(), _str.atoms.end(),
                       _str.k_vecs.begin(), _str.k_vecs.end(),
                       i_hold );

    Crystal::Structure::t_Atoms::iterator i_atom = _str.atoms.begin();
    Crystal::Structure::t_Atoms::iterator i_atom_end = _str.atoms.end();
    types :: t_int concx = 0, concy = 0; 
    i_hold = hold;
    for (; i_atom != i_atom_end; ++i_atom, ++i_hold)
    {
      if ( not ( i_atom->freeze & Crystal::Structure::t_Atom::FREEZE_T ) )
      {
        if ( std::abs( std::real(*i_hold) ) < types::tolerance )
          i_atom->type = rng.flip() ? 10.0*types::tolerance: -10.0*types::tolerance;
        else i_atom->type = std::real( *i_hold );
      }
      ( i_atom->type > 0.0 ) ? ++concx : --concx;

      ++i_atom;
      if ( not ( i_atom->freeze & Crystal::Structure::t_Atom::FREEZE_T ) )
      {
        if ( std::abs( std::imag(*i_hold) ) < types::tolerance )
          i_atom->type = rng.flip() ? 10.0*types::tolerance: -10.0*types::tolerance;
        else i_atom->type = std::imag( *i_hold );
      }
      ( i_atom->type > 0.0 ) ? ++concy : --concy;
    }

    // then normalize it while setting correct concentrations
    if ( single_c )
    {
      normalize( _str, 0, (types::t_real) concx - ( (types::t_real) N ) * x0 );
      normalize( _str, 1, (types::t_real) concy - ( (types::t_real) N ) * y0 );
      delete[] hold;
      return;
    }
    
    // Concentration is not set, but still constrained by load balancing
    // hence will randomly set concentration to ( x and load balanced y )
    // or ( y and load balanced x). The latter case may not always be possible. 
    x0 = (double) concx / (double) N;
    if ( rng.flip() or not can_inverse(x) )
    {  // x and load balanced y
      y0 = (double) concy / (double) N;
      x0 = get_x( y0 );
      normalize( _str, 1, 0 );
      normalize( _str, 0, (types::t_real) concx - ( (types::t_real) N ) * x0 );
      delete[] hold;
      return;
    }
     
    // y and load balanced x
    y0 = get_y( x0 );
    normalize( _str, 0, 0 );
    normalize( _str, 1, (types::t_real) concy - ( (types::t_real) N ) * y0 );
    delete[] hold;
  }

  void Concentration :: operator()( Object &_obj )
  {
    // computes concentrations first
    get( _obj );
    if( not single_c )
    {
      x0 = x; y0 = y;
      if ( rng.flip() or not can_inverse(x) ) x0 = get_x(y0);
      else                                    y0 = get_y(x0);
    }

    // compute number of changes to make per site type
    types::t_real to_change[2];
    to_change[0] = (types::t_real) N * ( x0 - x );
    to_change[1] = (types::t_real) N * ( y0 - y );
    if(     Fuzzy::gt(to_change[0],  -1.0)  and Fuzzy::le(to_change[0], 1.0 )
        and Fuzzy::gt(to_change[1],  -1.0)  and Fuzzy::le(to_change[1], 1.0 ) ) return;

    __ASSERT( sites.size() != _obj.Container().size(), "Inequivalent sizes.\n" )

    // Puts the positions which can be changed into a list
    std::vector<types::t_unsigned> pos[2];
    typedef Object :: t_Container :: const_iterator const_iterator;
    std::vector<bool> :: const_iterator i_site = sites.begin();
    const_iterator i_bit = _obj.Container().begin();
    const_iterator i_bit_end = _obj.Container().end();
    types::t_unsigned a = 0, b = 0;
    for(types::t_unsigned i=0; i_bit != i_bit_end; ++i_bit, ++i, ++i_site)
    {
      if( *i_site ) ++a;
      else ++b;
      types::t_unsigned site_index = *i_site ? 0: 1;
      if( to_change[site_index] > 0.0 and BitString::spin_down(*i_bit)  )
        pos[site_index].push_back( i );
      else if( to_change[site_index] < 0.0 and BitString::spin_up(*i_bit)  )
        pos[site_index].push_back( i );
    }

    // shuffle position lists
    types::t_unsigned (*ptr_to_rng)(const types::t_unsigned& )
        = &eo::random<types::t_unsigned>;
    for( types::t_unsigned i=0; i < 2; ++i )
    {
      if( Fuzzy::gt(to_change[i],  -1.0) and Fuzzy::le(to_change[i], 1.0 ) ) continue;
      std::vector<types::t_unsigned> :: iterator i_pos = pos[i].begin();
      std::vector<types::t_unsigned> :: iterator i_pos_end = pos[i].end();
      std::random_shuffle(i_pos, i_pos_end, ptr_to_rng );
      i_pos = pos[i].begin();
      for(; i_pos != i_pos_end; ++i_pos)
      {
        BitString::flip<Object::t_Container::value_type>(_obj.bitstring[*i_pos]);
        ( to_change[i] > 0 ) ? to_change[i] -= 2: to_change[i] += 2;
      
        // Fuzzy math at this point could create infinit loop
        if( Fuzzy::gt(to_change[i],  -1.0) and Fuzzy::leq(to_change[i], 1.0 ) ) break;
      }
      __DOASSERT( Fuzzy::leq(to_change[i],  -1.0) or Fuzzy::gt(to_change[i], 1.0 ),
                     "Concentration could not be set\n"
                  << "Incompatibility between required x/y and frozen atoms?\n"; )
    }

  // __DODEBUGCODE( get( _obj ); )
  // __DOASSERT( Fuzzy::neq(x, x0), x << " == " << x0 << "\n" )
  // __DOASSERT( Fuzzy::neq(y, y0), y << " == " << y0 << "\n" )
  }

  void Concentration :: normalize( Crystal::Structure &_str, const types::t_int _site, 
                                   types::t_real _tochange ) 
  {
    Crystal::Structure::t_Atoms::iterator i_end = _str.atoms.end();
    Crystal::Structure::t_Atoms::iterator i_which;
    Crystal::Structure::t_Atoms::iterator i_atom;
    while( _tochange < -1.0 or _tochange > 1.0 )
    {
      // first finds atom with type closest to zero from +1 or -1 side,
      // depending on _tochange
      i_atom = _str.atoms.begin();
      i_which = i_end;
      types::t_real minr = 0.0;
      for(; i_atom != i_end; ++i_atom )
      {
        if ( _site ) ++i_atom; 
        if ( i_atom->freeze & Crystal::Structure::t_Atom::FREEZE_T )
          goto endofloop;
        if ( _tochange > 0 )
        {
          if ( i_atom->type < 0 )
            goto endofloop;
          if ( minr != 0.0 and i_atom->type > minr )
            goto endofloop;
        }
        else // ( _tochange < 0 )
        {
          if ( i_atom->type > 0 )
            goto endofloop;
          if ( minr != 0.0 and i_atom->type < minr )
            goto endofloop;
        }

        i_which = i_atom;
        minr = i_atom->type;

endofloop: 
        if ( not _site ) ++i_atom;
      }
      __DOASSERT( i_which == i_end,
                  "Error while normalizing constituents of site " << _site << "\n")


      i_which->type = ( _tochange > 0 ) ? -1.0: 1.0;
      _tochange += ( _tochange > 0 ) ? -2: 2;
    }

    // finally, normalizes _str
    i_atom = _str.atoms.begin();
    for(; i_atom != i_end; ++i_atom )
    {
      if ( _site ) ++i_atom;
      i_atom->type = ( i_atom->type > 0 ) ? 1.0: -1.0;
      if ( not _site ) ++i_atom;
    }
    types::t_real concx = 0;
    types::t_real concy = 0;
    types :: t_unsigned N = _str.atoms.size() >> 1;
    i_atom = _str.atoms.begin();
    for(; i_atom != i_end; ++i_atom )
    {
      i_atom->type > 0 ? ++concx: --concx;
      ++i_atom;
      i_atom->type > 0 ? ++concy: --concy;
    }
    types::t_real result =  _site ? (types::t_real ) concy / (types::t_real) N:
                                    (types::t_real ) concx / (types::t_real) N;
    types::t_real inv = 2.0 / (types::t_real) N;
    __DOASSERT( std::abs( result - (_site ? y0:x0) ) > inv,
                   "Could not normalize site\n" << ( _site ? " x= ": " y= " )
                << result <<  ( _site ? " x0= ": " y0= " ) <<  ( _site ? x0: y0 )
                << "\n" )
  }

  void Concentration :: get( const Crystal::Structure &_str)
  {
    Crystal::Structure::t_Atoms :: const_iterator i_atom = _str.atoms.begin();
    Crystal::Structure::t_Atoms :: const_iterator i_atom_end = _str.atoms.end();
    types :: t_unsigned Nx = 0, Ny = 0;
    for(x = y = 0e0; i_atom != i_atom_end; ++i_atom )
      ( i_atom->site > 0 ) ? ++Ny, y += i_atom->type: 
                             ++Nx, x += i_atom->type;
    x /= (types::t_real) Nx;
    y /= (types::t_real) Ny;
  }

  void Concentration :: get( const Object &_obj )
  {
    __ASSERT( sites.size() != _obj.bitstring.size(),
              "Unequal sizes.\n" )

    // computes concentrations first
    Object::t_Container::const_iterator i_bit = _obj.bitstring.begin();
    Object::t_Container::const_iterator i_bit_end = _obj.bitstring.end();
    std::vector<bool>::const_iterator i_site = sites.begin();
    types::t_int concx = 0, concy = 0;
    for(; i_bit != i_bit_end; ++i_bit, ++i_site )
      if( *i_site ) *i_bit > 0 ? ++concx: --concx;
      else          *i_bit > 0 ? ++concy: --concy;

    // add frozen bits
    concx += Nfreeze_x;
    concy += Nfreeze_y;

    // finally normalizes
    x = (types::t_real) concx / (types::t_real) N;
    y = (types::t_real) concy / (types::t_real) N;
  }


  void Concentration :: setfrozen( const Crystal::Structure &_str )
  {
    __DOASSERT( _str.atoms.size() % 2 != 0,
                "Uneven number of sites in structure.\n" )
    N = _str.atoms.size() >> 1;

    Crystal::Structure::t_Atoms::const_iterator i_atom = _str.atoms.begin();
    Crystal::Structure::t_Atoms::const_iterator i_atom_end = _str.atoms.end();
    Nfreeze_x = 0; Nfreeze_y = 0;
    for(; i_atom != i_atom_end; ++i_atom )
    {
      if ( i_atom->freeze & Crystal::Structure::t_Atom::FREEZE_T )
        Nfreeze_x += i_atom->type > 0 ? 1 : -1; 
      else sites.push_back( true );
      ++i_atom;
      if ( i_atom->freeze & Crystal::Structure::t_Atom::FREEZE_T )
        Nfreeze_y += i_atom->type > 0 ? 1 : -1; 
      else sites.push_back( false );
    }

 //  if ( not _str.lattice ) return;
 //  if (    ( _str.lattice->sites[0].freeze & Crystal::Structure::t_Atom::FREEZE_T ) 
 //       or ( _str.lattice->sites[1].freeze & Crystal::Structure::t_Atom::FREEZE_T ) )
 //  {
 //    set_xy( x, y );
 //    Print::xmg << Print::Xmg::comment << " Setting Concentrations to x="
 //               << Print::fixed << Print::setprecision(3 ) << x 
 //               << " and y=" << y << Print::endl;
 //    if (    std::abs( x - get_x( y)) > types::tolerance 
 //         or std::abs( y - get_y( x)) > types::tolerance )
 //      Print::out << " WARNING: x and y pair are strain mismatched!! \n";
 //  }

  }


} // namespace pescan




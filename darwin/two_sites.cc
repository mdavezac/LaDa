//
//  Version: $Id$
//
#include <functional>
#include <algorithm>
#include <ext/algorithm>
#include <fstream>
#ifndef __PGI
  #include<ext/functional>
  using __gnu_cxx::compose1;
#else
  #include<functional>
  using std::compose1;
#endif

#include "two_sites.h"
#include "lamarck/atom.h"
#include "opt/compose_functors.h"

namespace TwoSites
{

  void rearrange_structure( Ising_CE::Structure &_str)
  {
    if ( not _str.lattice and _str.lattice->sites.size() != 2)
      return;

    std::list< Ising_CE::Structure::t_Atom > sites0;
    std::list< Ising_CE::Structure::t_Atom > sites1;
    Ising_CE::Structure::t_Atoms :: const_iterator i_atom = _str.atoms.begin();
    Ising_CE::Structure::t_Atoms :: const_iterator i_atom_end = _str.atoms.end();
    for(; i_atom != i_atom_end; ++i_atom )
      ( i_atom->site == 0 ) ?
        sites0.push_back( *i_atom ): sites1.push_back( *i_atom );

    std::list< Ising_CE::Structure::t_Atom > :: iterator i_0 = sites0.begin();
    std::list< Ising_CE::Structure::t_Atom > :: iterator i_end = sites0.end();
    std::list< Ising_CE::Structure::t_Atom > :: iterator i_1;
    atat::rVector3d translation = _str.lattice->sites[1].pos - _str.lattice->sites[0].pos; 
    types::t_real (*ptr_norm)(const atat::FixedVector<types::t_real, 3> &) = &atat::norm2;
    _str.atoms.clear();
    for(; i_0 != i_end; ++i_0 )
    {
      atat::rVector3d atom = i_0->pos + translation;
      i_1 = std::min_element ( sites1.begin(), sites1.end(), 
                opt::ref_compose2( std::less<types::t_real>(),
                                   compose1( std::ptr_fun(ptr_norm),
                                             bind2nd(std::minus<atat::rVector3d>(), atom) ),
                                     compose1( std::ptr_fun(ptr_norm),
                                               bind2nd(std::minus<atat::rVector3d>(), atom) ) ) );
      _str.atoms.push_back( *i_0 ); 
      _str.atoms.push_back( *i_1 ); 
    }

  }


  //  sets concentration from k-space values.
  //  individual need not be physical (ie S_i=+/-1) when fourier transformed to real-space
  //  a normalization procedure is applied which takes into account:
  //  (i) that this ga run is at set concentration (x =x0, y=y0)
  //  (ii) that x and y are only constrained by load-balancing
  void Concentration :: operator()( Ising_CE::Structure &_str )
  {
    types::t_complex  *hold = new types::t_complex[ N ];
    if ( not hold )
    {
      std::cerr << " Could not allocate memory in set_concentration!! " << std::endl; 
      exit(0);
    }

    // creates individual with unnormalized occupation
    types::t_complex  *i_hold = hold;
    TwoSites::Fourier( _str.atoms.begin(), _str.atoms.end(),
                       _str.k_vecs.begin(), _str.k_vecs.end(),
                       i_hold );

    Ising_CE::Structure::t_Atoms::iterator i_atom = _str.atoms.begin();
    Ising_CE::Structure::t_Atoms::iterator i_atom_end = _str.atoms.end();
    types :: t_int concx = 0, concy = 0; 
    i_hold = hold;
    for (; i_atom != i_atom_end; ++i_atom, ++i_hold)
    {
      if ( not ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) )
      {
        i_atom->type = std::real(*i_hold);
        if ( std::abs(i_atom->type) < types::tolerance )
          i_atom->type = rng.flip() ? 10.0*types::tolerance: -10.0*types::tolerance;
      }
      ( i_atom->type > 0 ) ? ++concx : --concx;

      ++i_atom;
      if ( not ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) )
      {
        i_atom->type = std::imag(*i_hold);
        if ( std::abs(i_atom->type) < types::tolerance )
          i_atom->type = rng.flip() ? 10.0*types::tolerance: -10.0*types::tolerance;
      }
      ( i_atom->type > 0 ) ? ++concy : --concy;
    }

    // then normalize it while setting correct concentrations
    if ( single_c )
    {
      normalize( _str, 0, (types::t_real) concx - ( (types::t_real) N ) * get_x(x) );
      normalize( _str, 1, (types::t_real) concy - ( (types::t_real) N ) * get_y(y) );
      delete[] hold;
      return;
    }
    
    // Concentration is not set, but still constrained by load balancing
    // hence will randomly set concentration to ( x and load balanced y )
    // or ( y and load balanced x). The latter case may not always be possible. 
    x = (double) concx / (double) N;
    if ( rng.flip() or not can_inverse(x) )
    {  // x and load balanced y
      y = (double) concy / (double) N;
      x = get_x( y );
      normalize( _str, 1, 0 );
      normalize( _str, 0, (types::t_real) concx - ( (types::t_real) N ) * x );
      delete[] hold;
      return;
    }
     
    // y and load balanced x
    y = get_y( x );
    normalize( _str, 0, 0 );
    normalize( _str, 1, (types::t_real) concy - ( (types::t_real) N ) * y );
    delete[] hold;
  }

  void Concentration :: operator()( Object &_obj )
  {
    // computes concentrations first
    set( _obj );
    if( not single_c )
    {
      x0 = x; y0 = y;
      if ( rng.flip() or not can_inverse(x) ) x0 = get_x(y0);
      else                                    y0 = get_y(x0);
    }

    // finally normalizes
    types::t_real xto_change = (types::t_real) N * ( x0 - x );
    types::t_real yto_change = (types::t_real) N * ( y0 - y );
    if (      xto_change > -1.0 and xto_change < 1.0 
         and  yto_change > -1.0 and xto_change < 1.0 ) return;
    do
    {
      types::t_unsigned i = rng.random(2*N-1);
      if ( sites[i] )
      {
        if ( xto_change > 1.0 and _obj.bitstring[i] < 0 )
          { _obj.bitstring[i] = 1; xto_change-=2; }
        else if ( xto_change < -1.0 and _obj.bitstring[i] > 0 )
          { _obj.bitstring[i] = -1; xto_change+=2; }
        continue;
      }
      
      if ( yto_change > 1.0 and _obj.bitstring[i] < 0 )
        { _obj.bitstring[i] = 1; yto_change-=2; }
      else if ( yto_change < -1.0 and _obj.bitstring[i] > 0 )
        { _obj.bitstring[i] = -1; yto_change+=2; }

    } while (    xto_change < -1.0 or xto_change > 1.0
              or yto_change < -1.0 or yto_change > 1.0 );
  }

  void Concentration :: normalize( Ising_CE::Structure &_str, const types::t_int _site, 
                                   types::t_real _tochange ) 
  {
    Ising_CE::Structure::t_Atoms::iterator i_end = _str.atoms.end();
    Ising_CE::Structure::t_Atoms::iterator i_which;
    Ising_CE::Structure::t_Atoms::iterator i_atom;
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
        if ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T )
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
      if ( i_which == i_end )
        throw std::runtime_error( "Error while normalizing x constituents\n" );

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
#ifdef _DEBUG
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
    if ( std::abs( result - (_site ? y:x) ) > inv )
      throw std::runtime_error("Could not normalize site\n" );
#endif
  }

  void Concentration :: set( const Ising_CE::Structure &_str)
  {
    Ising_CE::Structure::t_Atoms :: const_iterator i_atom = _str.atoms.begin();
    Ising_CE::Structure::t_Atoms :: const_iterator i_atom_end = _str.atoms.end();
    types :: t_unsigned Nx = 0, Ny = 0;
    for(; i_atom != i_atom_end; ++i_atom )
      ( i_atom->site > 0 ) ? ++Ny, y += i_atom->type: 
                             ++Nx, x += i_atom->type;
    x /= (types::t_real) Nx;
    y /= (types::t_real) Ny;
  }

  void Concentration :: set( const Object &_obj )
  {
    if ( not single_c ) return;
    if ( sites.size() != _obj.bitstring.size() ) return;

    // computes concentrations first
    Object::t_Container::const_iterator i_bit = _obj.bitstring.begin();
    Object::t_Container::const_iterator i_bit_end = _obj.bitstring.end();
    std::vector<bool>::const_iterator i_site = sites.begin();
    types::t_real concx = 0, concy = 0;
    for(; i_bit != i_bit_end; ++i_bit, ++i_site )
      if( *i_site ) concx += *i_bit > 0 ? 1: -1;
      else          concy += *i_bit > 0 ? 1: -1;

    // add frozen bits
    concx += Nfreeze_x;
    concy += Nfreeze_y;

    // finally normalizes
    x = (types::t_real) concx / (types::t_real) N;
    y = (types::t_real) concy / (types::t_real) N;
  }



} // namespace pescan




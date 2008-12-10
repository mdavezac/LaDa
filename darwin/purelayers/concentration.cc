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

#include "concentration.h"
#include "fourier.h"

namespace LaDa
{
  namespace GA
  {
    namespace Keepers
    {
      bool ConcTwo :: Load ( const TiXmlElement &_node )
      {
        double d;
      
        if ( not _node.Attribute("x", &d ) ) goto errorout;
        x = types::t_real(d);
        if ( not _node.Attribute("y", &d ) ) goto errorout;
        y = types::t_real(d);
      
        return true;
        errorout:
          std::cerr << "Could not Load BandGap::Keeper" << std::endl;
          return false;
      }
      bool ConcTwo :: Save( TiXmlElement &_node ) const
      {
        _node.SetDoubleAttribute("x", x );
        _node.SetDoubleAttribute("y", y );
      
        return true;
      }
    }

    namespace PureLayers
    {
  
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
        Fourier( _str.atoms.begin(), _str.atoms.end(),
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
      }
  
    } // namespace PureLayers


  } // namespace GA


} // namespace LaDa

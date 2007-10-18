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

#include <print/stdout.h>
#include <print/manip.h>
#include <lamarck/structure.h>
#include <lamarck/atom.h>
#include <opt/va_minimizer.h>

#include "molecularity.h"
#include "two_sites.h"

  // "Molecularity" means that we try and simulate growth conditions...
  // The input should be a 1x1xn supercell, with n either 100, 110 or 111 
  // each unit-cell is a "molecule", say InAs or GaSb, but not a mix of the two.


namespace Molecularity
{
  bool Evaluator :: Load( t_Individual &_indiv, const TiXmlElement &_node, bool _type )
  {
    t_Object &object = _indiv.Object();
    if ( not object.Load( _node ) ) return false;
    current_individual->quantities().front() =
       ( current_object->stress(0,0) + current_object->stress(1,1) ) * 0.5;
    current_individual->quantities().back() = current_object->cbm - current_object->vbm;

    return t_Base::Load( _indiv, _node, _type );
  }
  bool Evaluator :: Save( const t_Individual &_indiv, TiXmlElement &_node, bool _type ) const
  { 
    return _indiv.Object().Save(_node) and t_Base::Save( _indiv, _node, _type );
  }

  bool Evaluator :: Load( const TiXmlElement &_node )
  {
    if ( not lattice.Load( _node ) )
    {
      std::cerr << " Could not load lattice type from input!! " << std::endl; 
      return false;
    }
    Ising_CE::Structure::lattice = &lattice;
    if ( not structure.Load( _node ) )
    {
      std::cerr << " Could not load input structure!! " << std::endl; 
      return false;
    }
    if ( not structure.set_site_indices() )
    {
      std::cerr << " Could not set atomic indices! " << std::endl; 
      return false;
    }
    TwoSites::rearrange_structure(structure);
    if ( not consistency_check() )  return false;

    concentration.setfrozen( structure );


    if ( not vff.Load( _node ) )
    {
      std::cerr << " Could not load vff input!! " << std::endl; 
      return false;
    }
    if ( not pescan.Load( _node ) )
    {
      std::cerr << " Could not load pescan interface from input!! " << std::endl; 
      return false;
    }

    return true;
  }


  bool Evaluator :: consistency_check()
  {
    Ising_CE::Structure::t_Atoms :: iterator i_atom = structure.atoms.begin();
    Ising_CE::Structure::t_Atoms :: iterator i_atom_end = structure.atoms.end();
    for(; i_atom != i_atom_end; ++i_atom )
      if ( i_atom->site != 1 and i_atom->site != 0 ) return false;
    i_atom = structure.atoms.begin();

    for(; i_atom != i_atom_end; ++i_atom )
      if ( i_atom->site == 0 and not (i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) )
        break;
    if (     i_atom == i_atom_end
         and not (lattice.sites[0].freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) )
      lattice.sites[0].freeze |=  Ising_CE::Structure::t_Atom::FREEZE_T;
    if (     i_atom != i_atom_end 
         and (lattice.sites[0].freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) )
      for(i_atom = structure.atoms.begin(); i_atom != i_atom_end; i_atom+=2 )
        i_atom->freeze |= Ising_CE::Structure::t_Atom::FREEZE_T;
    if (  lattice.sites[0].freeze & Ising_CE::Structure::t_Atom::FREEZE_T )
    {
      std::cerr << "No atoms to optimize !? " << std::endl;
      return false;
    }

    // "freeze" type of second atom in unit-cell, so that it not placed within bitstring
    for(i_atom = structure.atoms.begin(); i_atom != i_atom_end; ++i_atom )
    {
      if ( i_atom->site != 1 ) continue;
      i_atom->freeze |= Ising_CE::Structure::t_Atom::FREEZE_T;
    }

    return true;
  }


  // Only need encode one atom per unit-cell
  void operator<<(Ising_CE::Structure &_str, const Object &_o)
  {
    Object::t_Container :: const_iterator i_var = _o.bitstring.begin();
    Object::t_Container :: const_iterator i_end = _o.bitstring.end();
    Ising_CE::Structure :: t_Atoms :: iterator i_atom = _str.atoms.begin();
    Ising_CE::Structure :: t_Atoms :: iterator i_atom_end = _str.atoms.end();
    for(; i_var != i_end and i_atom != i_atom_end; ++i_atom )
    {
      if ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) 
        { ++i_atom;  continue; }
      i_atom->type = *i_var > 0 ? 1.0 : -1.0;
      (++i_atom)->type = *i_var > 0 ? 1.0 : -1.0;
      ++i_var;
    }
  }
  void operator<<(Object &_o, const Ising_CE::Structure &_c)
  {
    _o.bitstring.clear(); _o.bitstring.reserve( _c.atoms.size() );
    Ising_CE::Structure :: t_Atoms :: const_iterator i_atom = _c.atoms.begin();
    Ising_CE::Structure :: t_Atoms :: const_iterator i_end = _c.atoms.end();
    for(; i_atom != i_end; i_atom += 2 )
      if ( not (i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T) )
        _o.bitstring.push_back( i_atom->type > 0 ? 1.0: -1.0 );
  }

  void Concentration :: operator()( Ising_CE::Structure &_str )
  {
    if ( not single_c ) return;

    types::t_complex  *hold = new types::t_complex[ N ];
    if ( not hold )
    {
      std::cerr << " Could not allocate memory in set_concentration!! " << std::endl; 
      exit(0);
    }

    // creates individual with unnormalized occupation
    types::t_complex  *i_hold = hold;
    Fourier( _str.atoms.begin(), _str.atoms.end(),
             _str.k_vecs.begin(), _str.k_vecs.end(),
            i_hold );

    Ising_CE::Structure::t_Atoms::iterator i_atom = _str.atoms.begin();
    Ising_CE::Structure::t_Atoms::iterator i_atom_end = _str.atoms.end();
    types :: t_int concx = 0;
    i_hold = hold;
    if( not single_c )
    {
      for (; i_atom != i_atom_end; i_atom += 2, ++i_hold)
        if ( not ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) )
        {
          if ( std::abs( std::real(*i_hold) ) < types::tolerance )
                i_atom->type = rng.flip() ? 1.0: -1.0;
          else  i_atom->type = std::real(*i_hold) > 0.0 ? 1.0: -1.0;
          (i_atom+1)->type = i_atom->type > 0.0 ? 1.0: -1.0;
        }
      return;
    }

    for (; i_atom != i_atom_end; i_atom+=2, ++i_hold)
    {
      if ( not ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) )
        i_atom->type = std::real(*i_hold);
      ( i_atom->type > 0 ) ? ++concx : --concx;
    }

    // then normalize it while setting correct concentrations
    normalize( _str, (types::t_real) concx - ( (types::t_real) N ) * x0 );
    delete[] hold;
    return;
  }

  // Takes an "unphysical" individual and set normalizes its sites _sites to +/-1,
  // after flipping the _tochange spins closest to zero.
  // ie sets the concentration
  void Concentration :: normalize( Ising_CE::Structure &_str, types::t_real _tochange ) 
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
      for(; i_atom != i_end; i_atom +=2 )
      {
        if ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) continue;
        if ( _tochange > 0 )
        {
          if ( i_atom->type < 0 )                    continue;
          if ( minr != 0.0 and i_atom->type > minr ) continue;
        }
        else // ( _tochange < 0 )
        {
          if ( i_atom->type > 0 )                    continue;
          if ( minr != 0.0 and i_atom->type < minr ) continue;
        }

        i_which = i_atom;
        minr = i_atom->type;
      }
      if ( i_which == i_end )
        throw std::runtime_error( "Error while normalizing x constituents\n" );

      i_which->type = ( _tochange > 0 ) ? -1.0: 1.0;
      _tochange += ( _tochange > 0 ) ? -2: 2;
    }

    // finally, normalizes _str
    i_atom = _str.atoms.begin();
    for(; i_atom != i_end; i_atom+=2 )
    {
      i_atom->type = ( i_atom->type > 0 ) ? 1.0: -1.0;
      (i_atom+1)->type = ( i_atom->type > 0 ) ? 1.0: -1.0;
    }
  }

  void Concentration :: set( const Ising_CE::Structure &_str)
  {
    Ising_CE::Structure::t_Atoms :: const_iterator i_atom = _str.atoms.begin();
    Ising_CE::Structure::t_Atoms :: const_iterator i_atom_end = _str.atoms.end();
    for(; i_atom != i_atom_end; i_atom +=2 )
      x += i_atom->type;
    x /= (types::t_real) N;
  }

  void Concentration :: setfrozen( const Ising_CE::Structure &_str )
  {
    // We need to redo this, since the number of "effective" atoms is haflf of what
    // SingleSite::Concentration expects.
    N = _str.atoms.size() >> 1;
    Ising_CE::Structure::t_Atoms::const_iterator i_atom = _str.atoms.begin();
    Ising_CE::Structure::t_Atoms::const_iterator i_atom_end = _str.atoms.end();
    Nfreeze = 0; 
    for(; i_atom != i_atom_end; i_atom+=2 )
      if ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T )
        Nfreeze += i_atom->type > 0 ? 1 : -1; 
  }
} // namespace Molecularity




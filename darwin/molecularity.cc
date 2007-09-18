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



  std::ostream& operator<<(std::ostream &_stream, const Object &_o)
  {
    Object::t_Container :: const_iterator i_var = _o.bitstring.begin();
    Object::t_Container :: const_iterator i_end = _o.bitstring.end();
    for(; i_var != i_end; ++i_var )
      _stream << ( *i_var > 0 ? '1' : '0' );
    _stream << " " << (Pescan::Keeper&) _o << " " 
            << ( _o.stress(0,0) + _o.stress(1,1) ) * 0.5  << " ";
    return _stream;
  }
  void operator<<(std::string &_str, const Object &_o)
  {
    std::ostringstream sstr;
    sstr << _o; _str = sstr.str();
  }
  void operator<<(Object &_o, const std::string &_c)
  {
    types::t_unsigned size = _c.size();
    _o.bitstring.resize( size );
    std::vector<types::t_real> :: iterator i_var = _o.bitstring.begin();
    std::vector<types::t_real> :: iterator i_end = _o.bitstring.end();
    for(types::t_unsigned n=0; i_var != i_end; ++i_var, ++n )
      *i_var = ( _c[n] == '1' ) ? 1.0: -1.0;
  }






  bool Evaluator :: Load( t_Individual &_indiv, const TiXmlElement &_node, bool _type )
  {
    t_Object &object = _indiv.Object();
    if ( not object.Load( _node ) ) return false;
    _indiv.quantities().back() = object.cbm - object.vbm; 
    _indiv.quantities().front() = ( object.stress(0,0) + object.stress(1,1) ) * 0.5;

    if ( _type == GA::LOADSAVE_SHORT )
    {
      if( not _node.Attribute("string") )
        return false;
      (_indiv.Object()) << std::string(_node.Attribute("string"));
      return true;
    }

    Ising_CE::Structure s; 
    if ( not s.Load(_node) )
      return false;
    (_indiv.Object()) << s;
    return true;
  }

  bool Evaluator :: Save( const t_Individual &_indiv, TiXmlElement &_node, bool _type ) const
  {
    if ( not _indiv.Object().Save( _node ) ) return false;
    if ( _type == GA::LOADSAVE_SHORT )
    {
      std::string str; str << _indiv.Object();
      _node.SetAttribute("string", str.c_str());
      return true;
    }

    Ising_CE::Structure s = structure; 
    s << _indiv.Object();
    Fourier( s.atoms.begin(),  s.atoms.end(),
             s.k_vecs.begin(), s.k_vecs.end() );
    s.print_xml(_node);

    return true;
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

    if ( not concentration.Load( _node ) ) 
    {
      std::cerr << " Could not load Concentration input!! " << std::endl; 
      return false;
    }
    concentration.N = structure.atoms.size() >> 1;
    Ising_CE::Structure::t_Atoms::const_iterator i_atom = structure.atoms.begin();
    Ising_CE::Structure::t_Atoms::const_iterator i_atom_end = structure.atoms.end();
    concentration.Nfreeze = 0; 
    for(; i_atom != i_atom_end; ++i_atom )
      if ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T )
        concentration.Nfreeze += i_atom->type > 0 ? 1 : -1; 
    
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
      if ( i_atom->site != 1 and i_atom->site != 0 )
        return false;
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


  bool Concentration :: Load ( const TiXmlElement &_node )
  {
    if ( not ( _node.Attribute("x") or _node.Attribute("y") ) )
    {
      single_c = false;
      return true;
    }
    if ( not _node.Attribute("x") ) return true;
    
    double d;
    _node.Attribute("x", &d);
    if( d < 0 or d > 1 ) goto errorout;
    single_c = true;
    x0 = 2.0 * (types::t_real) d - 1.0;
    return true;

errorout:
    std::cerr << "Error while reading concentration input\n";
    return false;
  }

  void Concentration :: operator()( Ising_CE::Structure &_str )
  {
    if ( not single_c ) return;

    types::t_unsigned N = (types::t_int) _str.atoms.size(); N = N>>1;
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

  void Concentration :: operator()( Object &_obj )
  {
    if ( not single_c ) return;

    // computes concentrations first
    Object::t_Container::const_iterator i_bit = _obj.bitstring.begin();
    Object::t_Container::const_iterator i_bit_end = _obj.bitstring.end();
    types::t_real conc = 0;
    for(; i_bit != i_bit_end; ++i_bit )
      conc += *i_bit > 0 ? 1: -1;

    // add frozen bits
    conc += Nfreeze;

    // finally normalizes
    x = conc / (types::t_real) N;
    types::t_real to_change = (types::t_real) N * x  - conc;
    if ( to_change > -1.0 and to_change < 1.0 ) return;
    do
    {
      types::t_unsigned i = rng.random( types::t_int(N-1) );
      if ( to_change > 1.0 and _obj.bitstring[i] < 0 )
        { _obj.bitstring[i] = 1; to_change-=2; }
      else if ( to_change < -1.0 and _obj.bitstring[i] > 0 )
        { _obj.bitstring[i] = -1; to_change+=2; }

    } while ( to_change < -1.0 or to_change > 1.0 );
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
    for(; i_atom != i_end; ++i_atom )
    {
      i_atom->type = ( i_atom->type > 0 ) ? 1.0: -1.0;
      (++i_atom)->type = ( i_atom->type > 0 ) ? 1.0: -1.0;
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
  void Concentration :: set( const Object &_obj )
  {
    if ( not single_c ) return;

    // computes concentrations first
    Object::t_Container::const_iterator i_bit = _obj.bitstring.begin();
    Object::t_Container::const_iterator i_bit_end = _obj.bitstring.end();
    types::t_real conc = 0;
    for(; i_bit != i_bit_end; ++i_bit )
      conc += *i_bit > 0 ? 1: -1;

    // add frozen bits
    conc += Nfreeze;

    // finally normalizes
    x = (types::t_real) conc / (types::t_real) N;
  }
} // namespace Molecularity




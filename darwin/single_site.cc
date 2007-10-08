//
//  Version: $Id$
//
#include <fstream>

#include <lamarck/atom.h>
#include <opt/va_minimizer.h>

#include "single_site.h"

namespace SingleSite
{
  std::ostream& operator<<(std::ostream &_stream, const Object &_o)
  {
    Object::t_Container :: const_iterator i_var = _o.bitstring.begin();
    Object::t_Container :: const_iterator i_end = _o.bitstring.end();
    for(; i_var != i_end; ++i_var )
      _stream << ( *i_var > 0 ? '1' : '0' );
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
  void operator<<(Ising_CE::Structure &_str, const Object &_o)
  {
    Object::t_Container :: const_iterator i_var = _o.bitstring.begin();
    Object::t_Container :: const_iterator i_end = _o.bitstring.end();
    Ising_CE::Structure :: t_Atoms :: iterator i_atom = _str.atoms.begin();
    Ising_CE::Structure :: t_Atoms :: iterator i_atom_end = _str.atoms.end();
    for(; i_var != i_end and i_atom != i_atom_end; ++i_atom )
    {
      if ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) 
        continue;
      i_atom->type = *i_var > 0 ? 1.0 : -1.0;
      ++i_var;
    }
  }
  void operator<<(Object &_o, const Ising_CE::Structure &_c)
  {
    _o.bitstring.clear(); _o.bitstring.reserve( _c.atoms.size() );
    Ising_CE::Structure :: t_Atoms :: const_iterator i_atom = _c.atoms.begin();
    Ising_CE::Structure :: t_Atoms :: const_iterator i_end = _c.atoms.end();
    for(; i_atom != i_end; ++i_atom )
      if ( not (i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T) )
        _o.bitstring.push_back( i_atom->type > 0 ? 1.0: -1.0 );
  }



  bool Concentration :: Load ( const TiXmlElement &_node )
  {
    std::string name = _node.Value();
    const TiXmlElement* parent = &_node;
    if ( name.compare("Concentration") )
      parent = _node.FirstChildElement("Concentration");
    if( not parent )
    {
      std::cerr << "Could not find Concentration Tag in input" << std::endl;
      return false;
    }
    const TiXmlAttribute *att = parent->FirstAttribute();
    for(; att ; att->Next() )
      LoadAttribute( *att );

    return true;
  }

  void Concentration :: LoadAttribute ( const TiXmlAttribute &_att )
  {
    single_c = false;
    std::string name = _att.Name();
    
    double d;
    d = _att.DoubleValue();
    if( d < 0 or d > 1 ) goto errorout;
    single_c = true;
    x0 = 2.0 * (types::t_real) d - 1.0;


errorout:
    std::cerr << "Error while reading concentration input\n";
  }


  void Concentration :: operator()( Ising_CE::Structure &_str )
  {
    types::t_unsigned N = (types::t_int) _str.atoms.size();
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
      for (; i_atom != i_atom_end; ++i_atom, ++i_hold)
        if ( not ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) )
        {
          if ( std::real( *i_hold ) < types::tolerance )
             i_atom->type = rng.flip() ? 1.0: -1.0;
          else 
           i_atom->type = std::real(*i_hold) > 0.0 ? 1.0: -1.0;
          std::cout << i_atom->type << " " << std::real( *i_hold ) << std::endl;
        }
      return;
    }

    for (; i_atom != i_atom_end; ++i_atom, ++i_hold)
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
    set( _obj );

    types::t_real to_change = (types::t_real) N * ( x0  - x );
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
      i_atom->type = ( i_atom->type > 0 ) ? 1.0: -1.0;
  }

  void Concentration :: set( const Ising_CE::Structure &_str)
  {
    Ising_CE::Structure::t_Atoms :: const_iterator i_atom = _str.atoms.begin();
    Ising_CE::Structure::t_Atoms :: const_iterator i_atom_end = _str.atoms.end();
    for(; i_atom != i_atom_end; ++i_atom )
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

  std::string Concentration :: print_out () const 
  {
    if ( not single_c ) return "Concentration Range";
    std::ostringstream sstr;
    sstr << "Single Concentration, x0 = " << x0;
    return sstr.str();
  }




} // SingleSite




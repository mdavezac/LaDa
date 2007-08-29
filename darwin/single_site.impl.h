//
//  Version: $Id$
//
#ifndef _SINGLE_SITE_IMPL_H_
#define _SINGLE_SITE_IMPL_H_

#include <algorithm>
#include <functional>
#include <ext/algorithm>
#ifndef __PGI
  #include<ext/functional>
  using __gnu_cxx::compose1;
#else
  #include<functional>
  using std::compose1;
#endif

#include "print/xmg.h"
#include "functors.h"
#include "concentration.h"

namespace SingleSite
{

  template<class T_INDIVIDUAL, class T_INDIV_TRAITS>
  bool Evaluator<T_INDIVIDUAL,T_INDIV_TRAITS> :: Load( t_Individual &_indiv,
                                                       const TiXmlElement &_node,
                                                       bool _type )
  {
    if ( _type == darwin::LOADSAVE_SHORT )
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

  template<class T_INDIVIDUAL, class T_INDIV_TRAITS>
  bool Evaluator<T_INDIVIDUAL,T_INDIV_TRAITS> :: Save( const t_Individual &_indiv,
                                                       TiXmlElement &_node,
                                                       bool _type ) const
  {
    if ( _type == darwin::LOADSAVE_SHORT )
    {
      std::string str; str << _indiv.Object();
      _node.SetAttribute("string", str.c_str());
      return true;
    }

    Ising_CE::Structure s = structure; 
    s << _indiv.Object();
    fourrier_to_kspace( s.atoms.begin(),  s.atoms.end(),
                         s.k_vecs.begin(), s.k_vecs.end() );
    s.print_xml(_node);

    return true;
  }

  //  sets concentration from k-space values.
  //  individual need not be physical (ie S_i=+/-1) when fourrier transformed to real-space
  //  a normalization procedure is applied which takes into account:
  //  (i) that this ga run is at set concentration (x =x0, y=y0)
  //  (ii) that x and y are only constrained by load-balancing
  template<class T_INDIVIDUAL, class T_INDIV_TRAITS>
  void Evaluator<T_INDIVIDUAL,T_INDIV_TRAITS> :: set_concentration( Ising_CE::Structure &_str )
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
    fourrier_to_rspace( _str.atoms.begin(), _str.atoms.end(),
                        _str.k_vecs.begin(), _str.k_vecs.end(),
                       i_hold );

    Ising_CE::Structure::t_Atoms::iterator i_atom = _str.atoms.begin();
    Ising_CE::Structure::t_Atoms::iterator i_atom_end = _str.atoms.end();
    types :: t_int concx = 0;
    i_hold = hold;
    for (; i_atom != i_atom_end; ++i_atom, ++i_hold)
    {
      if ( not ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) )
        i_atom->type = std::real(*i_hold);
      ( i_atom->type > 0 ) ? ++concx : --concx;
    }

    // then normalize it while setting correct concentrations
    normalize( _str, singlec ?   (types::t_real) concx - ( (types::t_real) N ) * x
                               : 0.0 );
    
    delete[] hold;
  }


  // Takes an "unphysical" individual and set normalizes its sites _sites to +/-1,
  // after flipping the _tochange spins closest to zero.
  // ie sets the concentration
  template<class T_INDIVIDUAL, class T_INDIV_TRAITS>
  void Evaluator<T_INDIVIDUAL,T_INDIV_TRAITS> :: normalize( Ising_CE::Structure &_str, 
                                                            types::t_real _tochange) 
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
        if ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T )
          continue;
        if ( _tochange > 0 )
        {
          if ( i_atom->type < 0 ) continue;
          if ( minr != 0.0 and i_atom->type > minr ) continue;
        }
        else // ( _tochange < 0 )
        {
          if ( i_atom->type > 0 ) continue;
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
#ifdef _DEBUG
    types::t_real concx = 0;
    types :: t_unsigned N = _str.atoms.size() >> 1;
    i_atom = _str.atoms.begin();
    for(; i_atom != i_end; ++i_atom )
      i_atom->type > 0 ? ++concx: --concx;
    types::t_real result =  (types::t_real ) concx / (types::t_real) N;
    types::t_real inv = 2.0 / (types::t_real) N;
    if ( std::abs( result - x) > inv )
      throw std::runtime_error("Could not normalize site\n" );
#endif
  }

  template<class T_INDIVIDUAL, class T_INDIV_TRAITS>
  bool Evaluator<T_INDIVIDUAL,T_INDIV_TRAITS> :: Load( const TiXmlElement &_node )
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
    singlec = false;
    const TiXmlElement *child = _node.FirstChildElement("Concentration");
    if ( not child ) return true;
    if ( not child->Attribute("x") ) return true;
    double d;
    child->Attribute("x", &d); x = (types::t_real) d;
    singlec = true;
    std::ostringstream sstr; 
    sstr << "Concentration is fixed to " << x;
    Print::xmg.add_comment( sstr.str() );

    return true;
  }

  template<class T_INDIVIDUAL, class T_INDIV_TRAITS>
  bool Evaluator<T_INDIVIDUAL,T_INDIV_TRAITS> :: Crossover ( t_Individual &_indiv1,
                                                             const t_Individual &_indiv2 )
  {
    t_Object &obj1 = _indiv1;
    const t_Object &obj2 = _indiv2;
    Object::t_Container :: iterator i_var1 = obj1.begin();
    Object::t_Container :: const_iterator i_var2 = obj2.begin();
    Object::t_Container :: const_iterator i_var2_end = obj2.end();
    for(; i_var2 != i_var2_end; ++i_var1, ++i_var2)
      if ( rng.flip(crossover_probability) ) 
        *i_var1 = *i_var2;
    structure << obj1;
    t_This::set_concentration( structure );
    obj1 << structure;
    return true;
  }

  // expects kspace value to exist!!
  template<class T_INDIVIDUAL, class T_INDIV_TRAITS>
  bool Evaluator<T_INDIVIDUAL,T_INDIV_TRAITS> :: Krossover( t_Individual  &_offspring,
                                                            const t_Individual &_parent,
                                                            bool _range )
  {
    t_Object &offspring  = _offspring;
    const t_Object &parent  = _parent;
    Ising_CE::Structure str1 = structure, str2 = structure;
    str1 << offspring; str2 << parent;
    fourrier_to_kspace( str1.atoms.begin(),  str1.atoms.end(),
                        str1.k_vecs.begin(), str1.k_vecs.end() );
    fourrier_to_kspace( str2.atoms.begin(),  str2.atoms.end(),
                        str2.k_vecs.begin(), str2.k_vecs.end() );
    if ( _range and str1.k_vecs.size() > 2 ) // range crossover ... kvec should be oredered according to size
    {  
      types::t_unsigned n = (types::t_unsigned)
         std::floor(   (types::t_real) rng.random ( str1.k_vecs.size() - 1 ) 
                     * (types::t_real) crossover_probability );
      __gnu_cxx::copy_n( str2.k_vecs.begin(), n, str1.k_vecs.begin() );
    }
    else // every point crossover
    {
      Ising_CE::Structure::t_kAtoms :: const_iterator i_p = str2.k_vecs.begin();
      Ising_CE::Structure::t_kAtoms :: const_iterator i_p_end = str2.k_vecs.end();
      Ising_CE::Structure::t_kAtoms :: iterator i_o = str1.k_vecs.begin();
      for ( ; i_p != i_p_end; ++i_p, ++i_o)
        if ( rng.flip(crossover_probability) ) 
          i_o->type = i_p->type;
    }
  
//   t_This::set_concentration( str1 );
    singlec ? ::set_concentration( str1, x ): ::set_concentration( str1 );
    offspring << str1;

    return true; // offspring has changed!
  }

  template<class T_INDIVIDUAL, class T_INDIV_TRAITS>
  eoGenOp<T_INDIVIDUAL>* Evaluator<T_INDIVIDUAL,T_INDIV_TRAITS> :: LoadGaOp(const TiXmlElement &_el )
  {
    std::string value = _el.Value();

    if ( value.compare("Crossover") == 0 )
    {
      _el.Attribute( "prob", &crossover_probability );
      crossover_probability = crossover_probability > 0 ? std::abs(crossover_probability) : 0.5;
      if ( crossover_probability > 1 ) crossover_probability = 0.5;
      std::ostringstream sstr;
      sstr << "Crossover rate = " << crossover_probability;
      Print::xmg.add_comment(sstr.str());
      // pointer is owned by caller !!
      return darwin::new_genop( *this, &t_This::Crossover, std::string( "Crossover" ) );
    }
    else if ( value.compare( "Krossover" ) == 0 )
    {
      bool att = false;
      _el.Attribute( "prob", &crossover_probability );
      crossover_probability = crossover_probability > 0 ? std::abs(crossover_probability) : 0.5;
      if ( crossover_probability > 1 ) crossover_probability = 0.5;
      std::ostringstream sstr;
      sstr << "Krossover, rate = " << crossover_probability;
      Print::xmg.add_comment(sstr.str());
      if ( _el.Attribute("type") )
      {
        std::string str =  _el.Attribute("type");
        if ( str.compare("range") == 0 ) 
          { att = true; Print::xmg.add_to_last( ", Range = true" ); }
      }
      // pointer is owned by caller !!
      return darwin::new_genop( *this, &t_This::Krossover, 
                                std::string( "Krossover" ), att);
    }

    return NULL;
  }

  template<class T_INDIVIDUAL, class T_INDIV_TRAITS>
  bool Evaluator<T_INDIVIDUAL,T_INDIV_TRAITS> :: Taboo(const t_Individual &_indiv )
  {
    if ( singlec ) return false;
    x = _indiv.get_concentration();
    return x > lessthan or x < morethan; // if true, _object is taboo
  }
  
  template<class T_INDIVIDUAL, class T_INDIV_TRAITS>
  darwin::Taboo_Base< T_INDIVIDUAL >* 
       Evaluator<T_INDIVIDUAL,T_INDIV_TRAITS> :: LoadTaboo(const TiXmlElement &_el )
  {
    if ( singlec ) return NULL;
    const TiXmlElement *child = _el.FirstChildElement( "Concentration" );
    if ( not child )
      return NULL;
    double d;
    if ( child->Attribute( "lessthan" ) )
      child->Attribute( "xlessthan", &d );
    lessthan = ( d > 0 and d < 1 ) ? 2.0*d-1.0: 1.0;
    if ( child->Attribute( "morethan" ) )
      child->Attribute( "morethan", &d );
    morethan = ( d > 0 and d < 1 ) ? 2.0*d-1.0: -1.0;
    if ( lessthan < morethan )
      return NULL;
   
    Print::xmg << Print::Xmg::comment << "Taboo x in [ " << 0.5*(morethan+1.0)
               << ", "  << 0.5*(lessthan+1.0) << "] " << Print::endl;
    // pointer is owned by caller !!
    return new darwin::TabooFunction< t_This >
                                    ( *this, &t_This::Taboo, "Taboo" );
  }

  template<class T_INDIVIDUAL, class T_INDIV_TRAITS>
  bool Evaluator<T_INDIVIDUAL,T_INDIV_TRAITS>::initialize( t_Individual &_indiv )
  {
    t_Object &object = _indiv.Object();
    object.bitstring.clear(); 
    Ising_CE::Structure::t_Atoms :: const_iterator i_atom = structure.atoms.begin();
    Ising_CE::Structure::t_Atoms :: const_iterator i_atom_end = structure.atoms.end();
    types::t_int concx = 0;
    for(; i_atom != i_atom_end; ++i_atom )
    {
      bool flip = rng.flip();
      if ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) 
        flip = ( i_atom->type > 0 );
      object.bitstring.push_back( flip ? -1.0: 1.0 );
      flip ? ++concx: --concx;
    }
    if ( singlec )
    {
      types::t_unsigned N = structure.atoms.size() >> 1; 
      types::t_real xto_change = (types::t_real) N * x  - concx;
      if ( xto_change > -1.0 and xto_change < 1.0 ) return true;
      do
      {
        types::t_unsigned i = rng.random(N-1);
        if ( xto_change > 1.0 and object.bitstring[i] < 0 )
          { object.bitstring[i] = 1; xto_change-=2; }
        else if ( xto_change < -1.0 and object.bitstring[i] > 0 )
          { object.bitstring[i] = -1; xto_change+=2; }

      } while ( xto_change < -1.0 or xto_change > 1.0 );
    }
    return true;
  }

} // namespace SingleSite
#endif // _TWOSITES_IMPL_H_

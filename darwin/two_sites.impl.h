//
//  Version: $Id$
//
#ifndef _TWOSITES_IMPL_H_
#define _TWOSITES_IMPL_H_

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
#include "print/stdout.h"
#include "functors.h"

namespace TwoSites
{
  template<class T_R_IT, class T_K_IT>
  void fourrier_to_kspace( T_R_IT _rfirst, T_R_IT _rend,
                           T_K_IT _kfirst, T_K_IT _kend ) // sets kvector values from rspace values
  {
    const std::complex<types::t_real>
       imath(0, -2*3.1415926535897932384626433832795028841971693993751058208);
    
    for (; _kfirst != _kend; ++_kfirst)
    {
      _kfirst->type = std::complex<types::t_real>(0);
      for(T_R_IT i_r( _rfirst ); i_r != _rend; ++i_r )
      {
        _kfirst->type +=   exp( imath * ( i_r->pos[0] * _kfirst->pos[0] +
                                          i_r->pos[1] * _kfirst->pos[1] +
                                          i_r->pos[2] * _kfirst->pos[2] ) )
                         * std::complex<types::t_real>(i_r->type, (++i_r)->type);
      }
    }
  }
  template<class T_R_IT, class T_K_IT, class T_O_IT >
  void fourrier_to_rspace( T_R_IT _rfirst, T_R_IT _rend,
                           T_K_IT _kfirst, T_K_IT _kend,
                           T_O_IT _rout ) // sets rvector values from kspace values
  {
    const std::complex<types::t_real>
       imath(0, 2*3.1415926535897932384626433832795028841971693993751058208);
    for (; _rfirst != _rend; _rfirst+=2, ++_rout)
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

  template<class T_INDIVIDUAL, class T_INDIV_TRAITS>
  bool Evaluator<T_INDIVIDUAL, T_INDIV_TRAITS> :: Load( t_Individual &_indiv,
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
  bool Evaluator<T_INDIVIDUAL, T_INDIV_TRAITS> :: Save( const t_Individual &_indiv, 
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
    TwoSites::fourrier_to_kspace( s.atoms.begin(),  s.atoms.end(),
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
    types::t_unsigned N = (types::t_int) _str.atoms.size(); N = N>>1;
    types::t_complex  *hold = new types::t_complex[ N ];
    if ( not hold )
    {
      std::cerr << " Could not allocate memory in set_concentration!! " << std::endl; 
      exit(0);
    }

    // creates individual with unnormalized occupation
    types::t_complex  *i_hold = hold;
    TwoSites::fourrier_to_rspace( _str.atoms.begin(), _str.atoms.end(),
                                  _str.k_vecs.begin(), _str.k_vecs.end(),
                                 i_hold );

    Ising_CE::Structure::t_Atoms::iterator i_atom = _str.atoms.begin();
    Ising_CE::Structure::t_Atoms::iterator i_atom_end = _str.atoms.end();
    types :: t_int concx = 0, concy = 0; 
    i_hold = hold;
    for (; i_atom != i_atom_end; ++i_atom, ++i_hold)
    {
      if ( not ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) )
        i_atom->type = std::real(*i_hold);
      ( i_atom->type > 0 ) ? ++concx : --concx;

      ++i_atom;
      if ( not ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) )
        i_atom->type = std::imag(*i_hold);
      ( i_atom->type > 0 ) ? ++concy : --concy;
    }

    // then normalize it while setting correct concentrations
    if ( x_vs_y.is_singlec() )
    {
      normalize( _str, 0, (types::t_real) concx - ( (types::t_real) N ) * x_vs_y.get_x(x) );
      normalize( _str, 1, (types::t_real) concy - ( (types::t_real) N ) * x_vs_y.get_y(y) );
      delete[] hold;
      return;
    }
    
    // Concentration is not set, but still constrained by load balancing
    // hence will randomly set concentration to ( x and load balanced y )
    // or ( y and load balanced x). The latter case may not always be possible. 
    x = (double) concx / (double) N;
    if ( rng.flip() or not x_vs_y.can_inverse(x) )
    {  // x and load balanced y
      y = (double) concy / (double) N;
      x = x_vs_y.get_x( y );
      normalize( _str, 1, 0 );
      normalize( _str, 0, (types::t_real) concx - ( (types::t_real) N ) * x );
      delete[] hold;
      return;
    }
     
    // y and load balanced x
    y = x_vs_y.get_y( x );
    normalize( _str, 0, 0 );
    normalize( _str, 1, (types::t_real) concy - ( (types::t_real) N ) * y );
    delete[] hold;
  }


  // Takes an "unphysical" individual and set normalizes its sites _sites to +/-1,
  // after flipping the _tochange spins closest to zero.
  // ie sets the concentration
  template<class T_INDIVIDUAL,class T_INDIV_TRAITS>
  void Evaluator<T_INDIVIDUAL,T_INDIV_TRAITS> :: normalize( Ising_CE::Structure &_str, 
                                                            const types::t_int _site, 
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
    if ( not structure.set_site_indices() )
    {
      std::cerr << " Could not set atomic indices! " << std::endl; 
      return false;
    }
    rearrange_structure(structure);
    if ( not consistency_check() )
      return false;
    if ( not x_vs_y.Load( _node ) )
    {
      std::cerr << " Could not load Concentration input!! " << std::endl; 
      return false;
    }
    if ( x_vs_y.is_singlec() )
      { x = x_vs_y.get_x(); y = x_vs_y.get_y(); }

    return true;
  }

  template<class T_INDIVIDUAL, class T_INDIV_TRAITS>
  bool Evaluator<T_INDIVIDUAL,T_INDIV_TRAITS> :: Crossover ( t_Individual &_indiv1,
                                                             const t_Individual &_indiv2 )
  {
    t_Object &obj1 = _indiv1;
    const t_Object &obj2 = _indiv2;
    Crossover<t_Object> :: crossover( crossover_rate );
    crossover( obj1, obj2 );
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
    TwoSites::fourrier_to_kspace( str1.atoms.begin(),  str1.atoms.end(),
                                  str1.k_vecs.begin(), str1.k_vecs.end() );
    TwoSites::fourrier_to_kspace( str2.atoms.begin(),  str2.atoms.end(),
                                  str2.k_vecs.begin(), str2.k_vecs.end() );

    // range crossover ... kvec should be oredered according to size
    if ( _range and str1.k_vecs.size() > 2 ) 
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
  
    t_This::set_concentration( str1 );
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
      Print::xmg << Print::Xmg::comment << "Crossover rate = "
                 << crossover_probability << Print::endl;
      // pointer is owned by caller !!
      return darwin::new_genop( *this, &t_This::Crossover, std::string( "Crossover" ) );
    }
    else if ( value.compare( "Krossover" ) == 0 )
    {
      bool att = false;
      _el.Attribute( "prob", &crossover_probability );
      crossover_probability = crossover_probability > 0 ? std::abs(crossover_probability) : 0.5;
      if ( crossover_probability > 1 ) crossover_probability = 0.5;
      Print::xmg << Print::Xmg::comment << "Krossover rate = " 
                 << crossover_probability << Print::endl;
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
    if ( x_vs_y.is_singlec() )
      return false;
    structure << (const t_Object&)_indiv;
    get_xy_concentrations( structure );
    return x > lessthan or x < morethan; // if true, _object is taboo
  }
  
  template<class T_INDIVIDUAL, class T_INDIV_TRAITS>
  darwin::Taboo_Base< T_INDIVIDUAL >* 
       Evaluator<T_INDIVIDUAL,T_INDIV_TRAITS> :: LoadTaboo(const TiXmlElement &_el )
  {
    if ( x_vs_y.is_singlec() )
      return NULL;
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
   
    Print::xmg << Print::Xmg::comment << Print::fixed << Print::setprecision(3) 
               << "Taboo x in [ " << 0.5*(morethan+1.0)
               << ", "  << 0.5*(lessthan+1.0) << "] " << Print::endl;
    // pointer is owned by caller !!
    return new darwin::TabooFunction< t_This >
                                    ( *this, &t_This::Taboo, "Taboo" );
  }

  template<class T_INDIVIDUAL, class T_INDIV_TRAITS>
  bool Evaluator<T_INDIVIDUAL,T_INDIV_TRAITS>::initialize( t_Object &_object )
  {
    _object.bitstring.clear(); 
    Ising_CE::Structure::t_Atoms :: const_iterator i_atom = structure.atoms.begin();
    Ising_CE::Structure::t_Atoms :: const_iterator i_atom_end = structure.atoms.end();
    types::t_int concx = 0;
    types::t_int concy = 0;
    for(; i_atom != i_atom_end; ++i_atom )
    {
      bool flip = rng.flip();
      if ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) 
        flip = ( i_atom->type > 0 );
      _object.bitstring.push_back( flip ? 1.0: -1.0 );
      flip ? ++concx: --concx;

      ++i_atom;
      flip = rng.flip();
      if ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) 
        flip = ( i_atom->type > 0 );
      _object.bitstring.push_back( flip ? 1.0: -1.0 );
      flip ? ++concy: --concy;
    }
    if ( x_vs_y.is_singlec() )
    {
      types::t_unsigned N = structure.atoms.size() >> 1; 
      types::t_real xto_change = (types::t_real) N * x  - concx;
      types::t_real yto_change = (types::t_real) N * y  - concy;
      if (      xto_change > -1.0 and xto_change < 1.0 
           and  yto_change > -1.0 and xto_change < 1.0 ) return true;
      do
      {
        types::t_unsigned i = 2 * rng.random(N-1);
        if ( xto_change > 1.0 and _object.bitstring[i] < 0 )
          { _object.bitstring[i] = 1; xto_change-=2; }
        else if ( xto_change < -1.0 and _object.bitstring[i] > 0 )
          { _object.bitstring[i] = -1; xto_change+=2; }
        
        if ( yto_change > -1.0 and yto_change < 1.0 ) continue;
        i = 2 * rng.random(N-1) + 1;
        if ( yto_change > 1.0 and _object.bitstring[i] < 0 )
          { _object.bitstring[i] = 1; yto_change-=2; }
        else if ( yto_change < -1.0 and _object.bitstring[i] > 0 )
          { _object.bitstring[i] = -1; yto_change+=2; }

      } while (    xto_change < -1.0 or xto_change > 1.0
                or yto_change < -1.0 or yto_change > 1.0 );
    }
    return true;
  }

  template<class T_INDIVIDUAL, class T_INDIV_TRAITS>
  void Evaluator<T_INDIVIDUAL,T_INDIV_TRAITS> :: get_xy_concentrations( const Ising_CE::Structure &_str)
  {
    Ising_CE::Structure::t_Atoms :: const_iterator i_atom = structure.atoms.begin();
    Ising_CE::Structure::t_Atoms :: const_iterator i_atom_end = structure.atoms.end();
    types :: t_unsigned Nx = 0, Ny = 0;
    for(; i_atom != i_atom_end; ++i_atom )
      ( i_atom->site > 0 ) ? ++Ny, y += i_atom->type: 
                             ++Nx, x += i_atom->type;
    x /= (types::t_real) Nx;
    y /= (types::t_real) Ny;
  }

  template<class T_INDIVIDUAL, class T_INDIV_TRAITS>
  bool Evaluator<T_INDIVIDUAL,T_INDIV_TRAITS> :: consistency_check()
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

    for(i_atom = structure.atoms.begin(); i_atom != i_atom_end; ++i_atom )
      if ( i_atom->site == 1 and not (i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) )
        break;
    if (     i_atom == i_atom_end
         and not (lattice.sites[1].freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) )
      lattice.sites[1].freeze |=  Ising_CE::Structure::t_Atom::FREEZE_T;
    if (     i_atom != i_atom_end
         and (lattice.sites[1].freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) )
      for(i_atom = structure.atoms.begin(); i_atom != i_atom_end; i_atom+=2 )
        (i_atom+1)->freeze |= Ising_CE::Structure::t_Atom::FREEZE_T;

    if (     ( lattice.sites[0].freeze & Ising_CE::Structure::t_Atom::FREEZE_T )
         and ( lattice.sites[1].freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) )
    {
      std::cerr << "No atoms to optimize !? " << std::endl;
      return false;
    }
    if (    ( lattice.sites[0].freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) 
         or ( lattice.sites[1].freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) )
    {
      get_xy_concentrations( structure );
      x_vs_y.set_xy( x, y );
      Print::xmg << Print::Xmg::comment << " Setting Concentrations to x="
                 << Print::fixed << Print::setprecision(3 ) << x 
                 << " and y=" << y << Print::endl;
      if (    std::abs(x - x_vs_y.get_x(y)) > types::tolerance 
           or std::abs(y - x_vs_y.get_y(x)) > types::tolerance )
        Print::out << " WARNING: x and y pair are strain mismatched!! \n";
    }

    return true;
  }


} // namespace TwoSites
#endif // _TWOSITES_IMPL_H_

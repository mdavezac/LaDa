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
    if ( not consistency_check() )  return false;

    if ( concentration.Load( _node ) ) return true;
    
    std::cerr << " Could not load Concentration input!! " << std::endl; 
    return false;
  }



  template<class T_INDIVIDUAL, class T_INDIV_TRAITS>
  eoGenOp<T_INDIVIDUAL>* Evaluator<T_INDIVIDUAL,T_INDIV_TRAITS> :: LoadGaOp(const TiXmlElement &_el )
  {
    std::string value = _el.Value();

    if ( value.compare("Crossover") == 0 )
    {
      Crossover<t_Individual, t_IndivTraits> *crossover
         = new Crossover<t_Individual, t_IndivTraits>(structure, concentration );
      if ( not crossover->Load( _el ) ) 
      {
        delete crossover;
        return NULL;
      }
      Print::xmg << Print::Xmg::comment << crossover->print_out() << Print::endl;
      // pointer is owned by caller !!
      return crossover;
    }
    else if ( value.compare("Mutation") == 0 )
    {
      Mutation<t_Individual, t_IndivTraits> *mutation
         = new Mutation<t_Individual, t_IndivTraits>(structure, concentration );
      if ( not mutation->Load( _el ) ) 
      {
        delete mutation;
        return NULL;
      }
      Print::xmg << Print::Xmg::comment << mutation->print_out() << Print::endl;
      // pointer is owned by caller !!
      return mutation;
    }
    else if ( value.compare( "Krossover" ) == 0 )
    {
      Krossover<t_Individual, t_IndivTraits> *krossover
         = new Krossover<t_Individual, t_IndivTraits>( concentration, structure );
      if ( not krossover->Load( _el ) ) 
      {
        delete krossover;
        return NULL;
      }
      Print::xmg << Print::Xmg::comment << krossover->print_out() << Print::endl;
      // pointer is owned by caller !!
      return krossover;
    }
    else if ( value.compare( "KMutation" ) == 0 )
    {
      KMutation<t_Individual, t_IndivTraits> *kmutation
         = new KMutation<t_Individual, t_IndivTraits>( concentration, structure );
      if ( not kmutation->Load( _el ) ) 
      {
        delete kmutation;
        return NULL;
      }
      Print::xmg << Print::Xmg::comment << kmutation->print_out() << Print::endl;
      // pointer is owned by caller !!
      return kmutation;
    }

    return NULL;
  }

  template<class T_INDIVIDUAL, class T_INDIV_TRAITS>
  bool Evaluator<T_INDIVIDUAL,T_INDIV_TRAITS> :: Taboo(const t_Individual &_indiv )
  {
    if ( concentration.is_singlec() )
      return false;
    structure << (const t_Object&)_indiv;
    concentration.set( structure );
    return concentration.x > lessthan or concentration.x < morethan; // if true, _object is taboo
  }
  
  template<class T_INDIVIDUAL, class T_INDIV_TRAITS>
  darwin::Taboo_Base< T_INDIVIDUAL >* 
       Evaluator<T_INDIVIDUAL,T_INDIV_TRAITS> :: LoadTaboo(const TiXmlElement &_el )
  {
    if ( concentration.is_singlec() )
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
  bool Evaluator<T_INDIVIDUAL,T_INDIV_TRAITS>::initialize( t_Individual &_indiv )
  {
    _indiv.Object().bitstring.clear(); 
    Ising_CE::Structure::t_Atoms :: const_iterator i_atom = structure.atoms.begin();
    Ising_CE::Structure::t_Atoms :: const_iterator i_atom_end = structure.atoms.end();
    types::t_int concx = 0;
    types::t_int concy = 0;
    for(; i_atom != i_atom_end; ++i_atom )
    {
      bool flip = rng.flip();
      if ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) 
        flip = ( i_atom->type > 0 );
      _indiv.Object().bitstring.push_back( flip ? 1.0: -1.0 );
      flip ? ++concx: --concx;

      ++i_atom;
      flip = rng.flip();
      if ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) 
        flip = ( i_atom->type > 0 );
      _indiv.Object().bitstring.push_back( flip ? 1.0: -1.0 );
      flip ? ++concy: --concy;
    }
    // sets bitstring to correct concentration if necessary
    concentration( structure, _indiv.Object(), concx, concy );
    return true;
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
      concentration.set_xy( concentration.x, concentration.y );
      Print::xmg << Print::Xmg::comment << " Setting Concentrations to x="
                 << Print::fixed << Print::setprecision(3 ) << concentration.x 
                 << " and y=" << concentration.y << Print::endl;
      if (    std::abs( concentration.x - concentration.get_x( concentration.y)) > types::tolerance 
           or std::abs( concentration.y - concentration.get_y( concentration.x)) > types::tolerance )
        Print::out << " WARNING: x and y pair are strain mismatched!! \n";
    }

    return true;
  }

  template<class T_INDIVIDUAL, class T_INDIV_TRAITS>
  bool Krossover<T_INDIVIDUAL,T_INDIV_TRAITS> :: Load( const TiXmlElement &_node )
  {
    std::string name = _node.Value();
    if ( name.compare("Krossover") ) return false;

    _node.Attribute( "prob", &rate );
    rate = rate > 0 ? std::abs(rate) : 0.5;
    if ( rate > 1 ) rate = 0.5;
    if ( _node.Attribute("type") )
    {
      std::string str =  _node.Attribute("type");
      if ( str.compare("range") == 0 )  do_range = true;
    }
    return true;
  }
  // expects kspace value to exist!!
  template<class T_INDIVIDUAL, class T_INDIV_TRAITS>
  void Krossover<T_INDIVIDUAL,T_INDIV_TRAITS> :: apply(eoPopulator<t_Individual> &_pop )
  {
    t_Object &offspring  = (*_pop).Object();
    const t_Object &parent  = _pop.select().Object();
    Ising_CE::Structure str1 = structure, str2 = structure;
    str1 << offspring; str2 << parent;
    TwoSites::fourrier_to_kspace( str1.atoms.begin(),  str1.atoms.end(),
                                  str1.k_vecs.begin(), str1.k_vecs.end() );
    TwoSites::fourrier_to_kspace( str2.atoms.begin(),  str2.atoms.end(),
                                  str2.k_vecs.begin(), str2.k_vecs.end() );

    // range crossover ... kvec should be oredered according to size
    if ( do_range and str1.k_vecs.size() > 2 ) 
    {  
      types::t_unsigned n = (types::t_unsigned)
         std::floor(   (types::t_real) rng.random ( str1.k_vecs.size() - 1 ) 
                     * (types::t_real) rate );
      __gnu_cxx::copy_n( str2.k_vecs.begin(), n, str1.k_vecs.begin() );
    }
    else // every point crossover
    {
      Ising_CE::Structure::t_kAtoms :: const_iterator i_p = str2.k_vecs.begin();
      Ising_CE::Structure::t_kAtoms :: const_iterator i_p_end = str2.k_vecs.end();
      Ising_CE::Structure::t_kAtoms :: iterator i_o = str1.k_vecs.begin();
      for ( ; i_p != i_p_end; ++i_p, ++i_o)
        if ( rng.flip(rate) ) 
          i_o->type = i_p->type;
    }
  
    concentration( str1 );
    offspring << str1;

    (*_pop).invalidate();
  }

  template<class T_INDIVIDUAL, class T_INDIV_TRAITS>
  bool KMutation<T_INDIVIDUAL,T_INDIV_TRAITS> :: Load( const TiXmlElement &_node )
  {
    std::string name = _node.Value();
    if ( name.compare("KMutation") ) return false;

    _node.Attribute( "prob", &rate );
    rate = rate > 0 ? std::abs(rate) : 0.5;
    if ( rate > 1 ) rate = 0.5;
    return true;
  }
  // expects kspace value to exist!!
  template<class T_INDIVIDUAL, class T_INDIV_TRAITS>
  void KMutation<T_INDIVIDUAL,T_INDIV_TRAITS> :: apply(eoPopulator<t_Individual> &_pop )
  {
    t_Object &object  = (*_pop).Object();
    Ising_CE::Structure str = structure;
    str << object;
    TwoSites::fourrier_to_kspace( str.atoms.begin(),  str.atoms.end(),
                                  str.k_vecs.begin(), str.k_vecs.end() );

    
    Ising_CE::Structure::t_kAtoms :: iterator i_k = str.k_vecs.begin();
    Ising_CE::Structure::t_kAtoms :: iterator i_k_end = str.k_vecs.end();
    types::t_real max = std::norm( i_k->type );
    for(; i_k != i_k_end; ++i_k)
      if ( max < std::norm( i_k->type ) ) max = std::norm( i_k->type );

    i_k = str.k_vecs.begin();
    max = std::sqrt( max );
    for(; i_k != i_k_end; ++i_k)
      if ( rng.flip(rate) )
        i_k->type = std::complex<types::t_real>( rng.uniform( max ), rng.uniform(max) );
  
    concentration( str );
    object << str;

    (*_pop).invalidate();
  }

} // namespace TwoSites
#endif // _TWOSITES_IMPL_H_

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


  template<class T_INDIVIDUAL>
  bool Evaluator<T_INDIVIDUAL> :: Load( t_Individual &_indiv,
                                        const TiXmlElement &_node,
                                        bool _type )
  {
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

  template<class T_INDIVIDUAL>
  bool Evaluator<T_INDIVIDUAL> :: Save( const t_Individual &_indiv,
                                        TiXmlElement &_node,
                                        bool _type ) const
  {
    if ( _type == GA::LOADSAVE_SHORT )
    {
      std::string str; str << _indiv.Object();
      _node.SetAttribute("string", str.c_str());
      return true;
    }

    Ising_CE::Structure s = structure; 
    s << _indiv.Object();
    t_FourierRtoK( s.atoms.begin(),  s.atoms.end(),
                   s.k_vecs.begin(), s.k_vecs.end() );
    s.print_xml(_node);

    return true;
  }

  template<class T_INDIVIDUAL>
  bool Evaluator<T_INDIVIDUAL> :: Load( const TiXmlElement &_node )
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
    
    concentration.setfrozen( structure );
    Print::xmg << Print::Xmg::comment << concentration.print() << Print::endl;

    return true;
  }

  template<class T_INDIVIDUAL>
  void Evaluator<T_INDIVIDUAL> :: presubmit( std::list<t_Individual> &_pop)
  {
    t_Individual pure;
    initialize( pure );
    // Pure A
    std::fill( pure.Object().Container().begin(), pure.Object().Container().end(), 1.0 );
    _pop.push_back(pure);
    // Pure B
    std::fill( pure.Object().Container().begin(), pure.Object().Container().end(), -1.0 );
    _pop.push_back(pure);
  } 
} // namespace SingleSite
#endif // _TWOSITES_IMPL_H_

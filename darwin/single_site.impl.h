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
  inline bool Evaluator<T_INDIVIDUAL> :: Load( const TiXmlElement &_node )
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


  template<class T_GAOPTRAITS, types::t_unsigned _D>
    typename Distance<T_GAOPTRAITS, _D>::t_ScalarFitnessQuantity 
      Distance<T_GAOPTRAITS, _D>::operator()( const t_Individual &_i1, 
                                              const t_Individual &_i2) const
      {
        typedef typename t_GATraits :: t_Quantity t_Quantity;
        typedef typename t_GATraits :: t_QuantityTraits t_QuantityTraits;
        if ( t_QuantityTraits::size( _i1.const_quantities() ) != 
             t_QuantityTraits::size( _i2.const_quantities() )     )
          throw std::runtime_error( "Individuals with differing number of quantities !?" );
        if ( t_QuantityTraits::size( _i1.const_quantities() ) < _D )
          throw std::runtime_error( "Inconsistent number of quantities !?" );
    
        concentration( _i1 ); 
        // Modifier::const_innermost() is a dirty hack which allows us to use
        // this same distance functor for all concentrations for which only x[0]
        // matters. 
        t_ScalarFitnessQuantity result = (t_ScalarFitnessQuantity)
          Modifier::const_innermost(concentration.x);
        concentration( _i2 ); 
        result -= (t_ScalarFitnessQuantity) Modifier::const_innermost(concentration.x);
        result = std::abs( result ) * xcoef;
    
        t_ScalarFitnessQuantity q(0);
        for( types::t_unsigned i = 0; i < _D; ++i )
          result +=   std::abs(  t_QuantityTraits::scalar( _i1.const_quantities(), i )
                               - t_QuantityTraits::scalar( _i1.const_quantities(), i ) )
                    * qcoefs[i];
        return result;
      }
  template<class T_GAOPTRAITS, types::t_unsigned _D>
    bool Distance<T_GAOPTRAITS, _D> :: Load( const TiXmlElement &_node )
    {
      bool result = false;
      double d = 0.0;
      if( _node.Attribute( "xcoef" ) ) _node.Attribute( "xcoef", &d );
      xcoef = (t_ScalarFitnessQuantity) d;
      if ( d ) result = true;
      for( types::t_unsigned i = 0; i < _D; ++i )
      {
        std::ostringstream sstr; 
        sstr << "qcoef" << i; 
        d = 0.0;
        if( _node.Attribute( sstr.str().c_str() ) )
          _node.Attribute( sstr.str().c_str(), &d );
        qcoefs[i] = d;
        if ( d ) result = true;
      } 
      return result;
    }

  template<class T_GATRAITS, types::t_int _D>
    Scaling::Base<T_GATRAITS>* new_Niche_from_xml( const TiXmlElement &_node )
    {
      if( not _node.Attribute("distance") ) return NULL;
      std::string name = _node.Attribute("distance");
      Scaling::Base<T_GATRAITS> *result;
      if( name != "Phenotypic" or name != "phenotypic" ) 
      {
        typedef Scaling::Niching<
            Scaling::Sharing::Triangular< Distance<T_GATRAITS, _D> > > t_Niche;
        t_Niche *result =  new t_Niche;
      }
      // Loading should be done by callee
      return result;
    }

} // namespace SingleSite
#endif // _TWOSITES_IMPL_H_

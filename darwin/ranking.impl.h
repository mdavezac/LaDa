//
//  Version: $Id$
//
#ifndef _RANKING_IMPL_H_
#define _RANKING_IMPL_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <tinyxml/tinyxml.h>

#include <opt/types.h>
#ifdef _MPI
#include <mpi/mpi_object.h>
#endif

namespace Ranking
{
  template<class T_GATRAITS>
  Container<T_GATRAITS> :: ~Container()
  {
    typename t_Container :: iterator i_ranking = container.begin();
    typename t_Container :: iterator i_end = container.end();
    for(; i_ranking != i_end; ++i_ranking )
    {
      delete (*i_ranking);
      *i_ranking = NULL;
    }
    container.clear();
  }

  template<class T_GATRAITS>
  void Container<T_GATRAITS> :: operator()(t_Population &_pop)
  {
    typename t_Container :: iterator i_ranking = container.begin();
    typename t_Container :: iterator i_end = container.end();
    for(; i_ranking != i_end; ++i_ranking )
      (*i_ranking)->operator()(_pop);
  }



  template<class T_SHARING>
  void Niching<T_SHARING> :: operator()(t_Population& _pop)
  {
    map( _pop );

    typename t_Population :: iterator i_indiv = _pop.begin();
    typename t_Population :: iterator i_end = _pop.end();
    typename t_Column :: const_iterator i_sum = sums.begin();
    for(; i_indiv != i_end; ++i_indiv );
      i_indiv->set_fitness( (t_ScalarFitnessQuantity&) i_indiv.fitness() / *i_sum );
  }

  template<class T_SHARING>
  void Niching<T_SHARING> :: map(const t_Population &_pop)
  {
    sums.clear(); sums.resize( _pop.size(), t_ScalarFitnessQuantity(0) );
    typename t_Population :: const_iterator i_2indiv;
    typename t_Population :: const_iterator i_indiv = _pop.begin();
    typename t_Population :: const_iterator i_end = _pop.end();
    typename t_Column :: iterator i_sum = sums.begin();
    typename t_Column :: iterator i_2sum;
    for(; i_indiv != i_end; ++i_indiv, ++i_sum );
    {
      i_2indiv = i_indiv + 1;
      i_2sum = i_sum + 1;
      *i_sum += t_ScalarFitnessQuantity(1);
      for(; i_2indiv != i_end; ++i_2indiv, ++i_2sum );
      {
        t_ScalarFitnessQuantity d = share(*i_indiv, *i_2indiv);
        (*i_sum)  += d; (*i_2sum) += d;
      }
    }
  }

  template<class T_SHARING>
  void ParetoRanking<T_SHARING> :: operator()(t_Population& _pop)
  {
    map( _pop );

    typename t_Population :: iterator i_indiv = _pop.begin();
    typename t_Population :: iterator i_end = _pop.end();
    typename t_Column :: const_iterator i_sum = sums.begin();
    for(; i_indiv != i_end; ++i_indiv );
      i_indiv->set_fitness( (t_ScalarFitnessQuantity)*i_sum );
  }

  template<class T_SHARING>
  void ParetoRanking<T_SHARING> :: map(const t_Population &_pop)
  {
    sums.clear(); sums.resize( _pop.size(), 0 );
    typename t_Population :: const_iterator i_2indiv;
    typename t_Population :: const_iterator i_indiv = _pop.begin();
    typename t_Population :: const_iterator i_end = _pop.end();
    t_Column :: iterator i_sum = sums.begin();
    t_Column :: iterator i_2sum;
    for(; i_indiv != i_end; ++i_indiv, ++i_sum );
    {
      i_2indiv = i_indiv + 1;  
      i_2sum = i_sum + 1;
      for(; i_2indiv != i_end; ++i_2indiv, ++i_2sum );
      {
        // minimizes by default -- see objectives.h
        if( i_indiv->fitness()  > i_2indiv->fitness() ) ++(*i_sum);
        if( i_2indiv->fitness() > i_indiv->fitness()  ) ++(*i_2sum);
      }
    }
  }



} // namespace Ranking

namespace Distance
{
  template<class T_GAOPTRAITS>
  typename Hamming<T_GAOPTRAITS>::t_ScalarFitnessQuantity 
     Hamming<T_GAOPTRAITS>::operator()( const t_Individual &_i1, const t_Individual &_i2) const
     {
       if ( _i1.Object().Container().size() != _i1.Object().Container().size() )
         throw std::runtime_error("individuals of different size in Distance::Hamming!\n");
       typename t_Object::const_iterator i_bit2 = _i2.Object().begin();
       typename t_Object::const_iterator i_bit1 = _i1.Object().begin();
       typename t_Object::const_iterator i_bit1_end = _i1.Object().end();
       t_ScalarFitnessQuantity result;
       for(; i_bit1 != i_bit1_end; ++i_bit1, ++i_bit2 )
         result += std::abs(   t_ScalarFitnessQuantity(*i_bit1) 
                             - t_ScalarFitnessQuantity(*i_bit2) );
       return result;
     }
} // namespace Distance

namespace Sharing
{
  template <class T_DISTANCE>
  bool Triangular<T_DISTANCE> :: Load( const TiXmlElement &_node )
  {
    // First finds correct node
    std::string name = _node.Value();
    const TiXmlElement *parent;
    if( name == "Sharing" )
      parent = &_node;
    else
      parent = _node.FirstChildElement("Sharing");
    
    for(; parent; parent = parent->NextSiblingElement("Sharing") )
    {
      if ( not parent->Attribute("type") ) continue;
      name = parent->Attribute("type");
      if(    name == "Triangular"
          or name == "triangular" ) break;
    }

    if ( not parent ) return false;
    if ( parent->Attribute("alpha") )
    {
      int d = 0;
      parent->Attribute("alpha", &d );
      if ( d < 0 ) 
      {
        std::cerr << "Invalid alpha = " << d 
                  << " in Sharing::Triangular " << std::endl;
        return false;
      } 
        
      alpha = types::t_unsigned( std::abs(d) );
    }
    if ( parent->Attribute("sigma") )
    {
      double d = 0;
      parent->Attribute("sigma", &d );
      if ( d < 0 ) 
      {
        std::cerr << "Invalid sigma = " << d 
                  << " in Sharing::Triangular " << std::endl;
        return false;
      } 
        
      sigma = t_ScalarFitnessQuantity( d );
    }
  }
  
  template<class T_DISTANCE>
  std::string Triangular<T_DISTANCE> :: print_out() const
  { 
    std::ostringstream sstr; 
    sstr << "Triangular Sharing, alpha=" << alpha << ", sigma=" << sigma
         << ", distance: " << distance.print_out();
    return sstr.str();
  }

  template<class T_DISTANCE>
    typename Triangular<T_DISTANCE>::t_ScalarFitnessQuantity 
       Triangular<T_DISTANCE> :: operator()(const t_Individual& _i1, 
                                            const t_Individual& _i2 ) const
    {
      t_ScalarFitnessQuantity d = distance( _i1, _i2 );
      if ( d > sigma ) return t_ScalarFitnessQuantity(0);
      if( alpha == -1 )
        return t_ScalarFitnessQuantity(1) - d / sigma;
      return std::pow( t_ScalarFitnessQuantity(1) - d / sigma, alpha );
    }
                            
} // namespace Sharing


namespace Ranking
{
  template<class T_GATRAITS> 
  Base<T_GATRAITS>* new_from_xml( const TiXmlElement &_node )
  {
    const TiXmlElement *parent = &_node;
    std::string name = parent->Value();
    if( name != "Ranking" )
      parent = _node.FirstChildElement("Ranking");
    if( not parent ) return NULL;

    Container<T_GATRAITS> *container = new Container<T_GATRAITS>;

    for(; parent; parent = parent->NextSiblingElement("Ranking") )
    {
      if( not parent->Attribute("type") ) return false;
      name = parent->Attribute("type");

      if( name == "Niching" or name == "niching" )
        { container.push_back( new_Niche_from_xml<T_GATRAITS>( *parent ) ); }
    }

    Base<T_GATRAITS>* result = NULL;
    switch( container->size() )
    {
      case 1: result = container->front(); 
      case 0: delete container; break;
      default: result = dynamic_cast< Base<T_GATRAITS>* >( container ); break;
    }

    return result;
  }

  template<class T_GATRAITS> 
  Base<T_GATRAITS>* new_Niche_from_xml( const TiXmlElement &_node )
  {
    typedef Niching< Sharing::Triangular< Distance::Hamming<T_GATRAITS> > > t_Niche;
    t_Niche result =  new t_Niche;
    if ( result and result->Load(_node ) ) return result;
    if ( result ) delete result;
    return NULL;

//   if( not _node.Attribute("sharing") ) return false;
//   name = _node.Attribute("sharing");
//   if( name != "Triangular" and name != "triangular" ) return false;
//
//   if( not _node.Attribute("distance") ) return false;
//   name = _node.Attribute("distance");
//   if( name != "Hamming" and name != "hamming" ) return false;
  }
}

#endif // _RANKING_H_

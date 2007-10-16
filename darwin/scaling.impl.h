//
//  Version: $Id$
//
#ifndef _RANKING_IMPL_H_
#define _RANKING_IMPL_H_

namespace Scaling
{
  template<class T_GATRAITS>
  Container<T_GATRAITS> :: ~Container()
  {
    typename t_Container :: iterator i_ranking = rankers.begin();
    typename t_Container :: iterator i_end = rankers.end();
    for(; i_ranking != i_end; ++i_ranking )
    {
      delete (*i_ranking);
      *i_ranking = NULL;
    }
    rankers.clear();
  }

  template<class T_GATRAITS>
  inline std::string Container<T_GATRAITS> :: what_is() const
  {
    std::ostringstream sstr;
    sstr << "Scaling begin{ ";
    typename t_Container::const_iterator i_rankers = rankers.begin();
    typename t_Container::const_iterator i_end = rankers.begin();
    for(; i_rankers != i_end; ++i_rankers)
      sstr << (*i_rankers)->what_is() << " "; 
    sstr << "} end"; 
    return  sstr.str();
  }

  template<class T_GATRAITS>
  inline void Container<T_GATRAITS> :: operator()(t_Population &_pop)
  {
    typename t_Container :: iterator i_ranking = rankers.begin();
    typename t_Container :: iterator i_end = rankers.end();
    for(; i_ranking != i_end; ++i_ranking )
      (*i_ranking)->operator()(_pop);
  }

  template<class T_GATRAITS>
   inline Base<T_GATRAITS>* Container<T_GATRAITS> :: pop_front()
   { 
     Base<T_GATRAITS> *copy = rankers.front();
     rankers.pop_front();
     return copy;
   }


  template<class T_SHARING>
  void Niching<T_SHARING> :: operator()(t_Population& _pop)
  {
    map( _pop );

    typename t_Population :: iterator i_indiv = _pop.begin();
    typename t_Population :: iterator i_end = _pop.end();
    typename t_Column :: const_iterator i_sum = sums.begin();
    for(; i_indiv != i_end; ++i_indiv );
    {
      t_ScalarFitness& fitness = i_indiv->fitness();
      fitness = ( (t_ScalarFitnessQuantity&) fitness ) / *i_sum;
    }
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
        t_ScalarFitnessQuantity d = sharing(*i_indiv, *i_2indiv);
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
      --(*i_sum);
      for(; i_2indiv != i_end; ++i_2indiv, ++i_2sum );
      {
        // minimizes by default -- see objectives.h
        if     ( i_indiv->fitness()   >  i_2indiv->fitness() )   --(*i_sum);
        else if( i_2indiv->fitness()  > i_indiv->fitness()   )               --(*i_2sum);
        else if( i_indiv->fitness()  == i_2indiv->fitness()  ) { --(*i_sum); --(*i_2sum); }
      }
    }
  }


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
      if ( parent->Attribute("d0") )
      {
        double d = 0;
        parent->Attribute("d0", &d );
        if ( d < 0 ) 
        {
          std::cerr << "Invalid d0 = " << d 
                    << " in Sharing::Triangular " << std::endl;
          return false;
        } 
          
        d_0 = t_ScalarFitnessQuantity( d );
      }
      return true;
    }
    
    template<class T_DISTANCE>
    std::string Triangular<T_DISTANCE> :: what_is() const
    { 
      std::ostringstream sstr; 
      sstr << "Sharing::Triangular( alpha=" << alpha << ", d0=" << d_0
           << ", distance: " << distance.what_is() << " )";
      return sstr.str();
    }

    template<class T_DISTANCE>
      typename Triangular<T_DISTANCE>::t_ScalarFitnessQuantity 
         Triangular<T_DISTANCE> :: operator()(const t_Individual& _i1, 
                                              const t_Individual& _i2 ) const
      {
        t_ScalarFitnessQuantity d = distance( _i1, _i2 );
        if ( d > d_0 ) return t_ScalarFitnessQuantity(0);
        if( alpha == 1 ) return t_ScalarFitnessQuantity(1) - d / d_0;
        return std::pow( t_ScalarFitnessQuantity(1) - d / d_0, (double) alpha );
      }
                              
  } // namespace Sharing



  namespace Distance
  {
    template<class T_GAOPTRAITS>
    typename GeneralHamming<T_GAOPTRAITS>::t_ScalarFitnessQuantity 
       GeneralHamming<T_GAOPTRAITS>::operator()( const t_Individual &_i1, const t_Individual &_i2) const
       {
         if ( _i1.Object().Container().size() != _i1.Object().Container().size() )
           throw std::runtime_error("individuals of different size in Distance::GeneralHamming!\n");
         typename t_Object::const_iterator i_bit2 = _i2.Object().begin();
         typename t_Object::const_iterator i_bit1 = _i1.Object().begin();
         typename t_Object::const_iterator i_bit1_end = _i1.Object().end();
         t_ScalarFitnessQuantity result;
         for(; i_bit1 != i_bit1_end; ++i_bit1, ++i_bit2 )
           result += std::abs(   t_ScalarFitnessQuantity(*i_bit1) 
                               - t_ScalarFitnessQuantity(*i_bit2) );
         return result;
       }

    template<class T_GAOPTRAITS>
    typename Hamming<T_GAOPTRAITS>::t_ScalarFitnessQuantity 
       Hamming<T_GAOPTRAITS>::operator()( const t_Individual &_i1, const t_Individual &_i2) const
       {
         typedef typename t_GATraits::t_QuantityTraits::t_ScalarQuantityTraits t_SQTraits;
         if ( _i1.Object().Container().size() != _i1.Object().Container().size() )
           throw std::runtime_error("individuals of different size in Distance::Hamming!\n");
         typename t_Object::const_iterator i_bit2 = _i2.Object().begin();
         typename t_Object::const_iterator i_bit1 = _i1.Object().begin();
         typename t_Object::const_iterator i_bit1_end = _i1.Object().end();
         t_ScalarFitnessQuantity result;
         for(; i_bit1 != i_bit1_end; ++i_bit1, ++i_bit2 )
           result +=  t_ScalarFitnessQuantity(t_SQTraits::equal( *i_bit1, *i_bit2) ? 1: 0);
         return result;
       }
  } // namespace Distance


  template<class T_GATRAITS> 
  Base<T_GATRAITS>* new_from_xml( const TiXmlElement &_node )
  {
    const TiXmlElement *parent = &_node;
    std::string name = parent->Value();
    if( name != "Scaling" )
      parent = _node.FirstChildElement("Scaling");
    if( not parent ) return NULL;

    Container<T_GATRAITS> *container = new Container<T_GATRAITS>;

    for(; parent; parent = parent->NextSiblingElement("Scaling") )
    {
      if( not parent->Attribute("type") ) return false;
      name = parent->Attribute("type");

      if( name == "Niching" or name == "niching" )
        { container->push_back( new_Niche_from_xml<T_GATRAITS>( *parent ) ); }
      if( T_GATRAITS::t_QuantityTraits::is_vector and ( name == "Pareto" or name == "pareto" ) )
        { container->push_back( new ParetoRanking<T_GATRAITS>() ); }
    }

    Base<T_GATRAITS>* result = NULL;
    switch( container->size() )
    {
      case 1: result = container->pop_front(); 
      case 0: delete container; break;
      default: result = dynamic_cast< Base<T_GATRAITS>* >( container ); break;
    }

    return result;
  }

  template<class T_GATRAITS> 
  Base<T_GATRAITS>* new_Niche_from_xml( const TiXmlElement &_node )
  {
    if( not _node.Attribute("distance") ) return NULL;
    std::string name = _node.Attribute("distance");
    if( name != "GeneralHamming" and name != "generalhamming" ) 
    {
      typedef Niching< Sharing::Triangular< Distance::GeneralHamming<T_GATRAITS> > > t_Niche;
      t_Niche *result =  new t_Niche;
      if ( result and result->Load(_node ) ) return result;
      if ( result ) delete result;
      return NULL;
    }

    if( name != "Hamming" and name != "hamming" ) return NULL;
    typedef Niching< Sharing::Triangular< Distance::Hamming<T_GATRAITS> > > t_Niche;
    t_Niche *result =  new t_Niche;
    if ( result and result->Load(_node ) ) return result;
    if ( result ) delete result;
    return NULL;
  }

} // namespace Scaling
#endif // _RANKING_IMPL_H_

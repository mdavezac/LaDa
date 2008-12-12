//
//  Version: $Id$
//
#ifndef _RANKING_IMPL_H_
#define _RANKING_IMPL_H_

#include <boost/regex.hpp>

#include <opt/fuzzy.h>

namespace LaDa
{
  namespace Scaling
  {

    template<class T_GATRAITS> 
    inline void Base<T_GATRAITS> :: operator()( typename T_GATRAITS::t_Population& _1,
                                                typename T_GATRAITS::t_Population& _2 )
    {
      GA::Aggregate< typename T_GATRAITS::t_Population > aggregate( _1, _2); 
      operator()( aggregate ); 
    }


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
      typename t_Container::const_iterator i_end = rankers.end();
      for(; i_rankers != i_end; ++i_rankers)
        sstr << (*i_rankers)->what_is() << " "; 
      sstr << " }end"; 
      return  sstr.str();
    }

    template<class T_GATRAITS> template<class TPOPULATION >
    inline void Container<T_GATRAITS> :: _operator_( TPOPULATION &_pop)
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


    template<class T_SHARING> template<class TPOPULATION >
    inline void Niching<T_SHARING> :: _operator_( TPOPULATION &_pop)
    {
      map( _pop );

      typename TPOPULATION :: iterator i_indiv = _pop.begin();
      typename TPOPULATION :: iterator i_end = _pop.end();
      typename t_Column :: const_iterator i_sum = sums.begin();
      for(; i_indiv != i_end; ++i_indiv, ++i_sum )
      {
        t_Individual& indiv  = Modifier::innermost( *i_indiv );
        t_ScalarFitnessQuantity fitness = indiv.fitness();
        indiv.set_fitness( ( fitness - offset ) / *i_sum );
      }
    }

    template<class T_SHARING> template< class TPOPULATION >
    void Niching<T_SHARING> :: map(const TPOPULATION &_pop)
    {
      sums.clear(); sums.resize( _pop.size(), t_ScalarFitnessQuantity(0) );
      typename TPOPULATION :: const_iterator i_2indiv;
      typename TPOPULATION :: const_iterator i_indiv = _pop.begin();
      typename TPOPULATION :: const_iterator i_end = _pop.end();
      typename t_Column :: iterator i_sum = sums.begin();
      typename t_Column :: iterator i_2sum;
      offset = t_ScalarFitnessQuantity(0);
      for(; i_indiv != i_end; ++i_indiv, ++i_sum )
      {
        const t_Individual& indiv  = Modifier::const_innermost( *i_indiv );
        // Finds the offset by which all fitnesses should be shifted
        // to make them negative
        if ( (const t_ScalarFitnessQuantity& ) indiv.fitness() > offset) 
          offset = (const t_ScalarFitnessQuantity& ) indiv.fitness();

        i_2indiv = i_indiv + 1;
        i_2sum = i_sum + 1;
        *i_sum += t_ScalarFitnessQuantity(1);
        for(; i_2indiv != i_end; ++i_2indiv, ++i_2sum )
        {
          t_ScalarFitnessQuantity d;
          d = sharing( indiv, Modifier::const_innermost(*i_2indiv) );
          (*i_sum)  += d; (*i_2sum) += d;
        }
      }

      // if offset is negative, then no offset is needed.
      if ( offset < t_ScalarFitnessQuantity(0) ) 
        offset = t_ScalarFitnessQuantity(0);
    }

    template<class T_SHARING> template< class TPOPULATION >
    inline void ParetoRanking<T_SHARING> :: _operator_( TPOPULATION & _pop)
    {
      map( _pop );

      typename TPOPULATION :: iterator i_indiv = _pop.begin();
      typename TPOPULATION :: iterator i_end = _pop.end();
      t_Column :: const_iterator i_sum = sums.begin();
      for(; i_indiv != i_end; ++i_indiv, ++i_sum )
        Modifier::innermost(*i_indiv).set_fitness(*i_sum);
    }

    template<class T_SHARING> template< class TPOPULATION >
    void ParetoRanking<T_SHARING> :: map(const TPOPULATION &_pop)
    {
      sums.clear(); sums.resize( _pop.size(), 0 );
      typename TPOPULATION :: const_iterator i_2indiv;
      typename TPOPULATION :: const_iterator i_indiv = _pop.begin();
      typename TPOPULATION :: const_iterator i_end = _pop.end();
      t_Column :: iterator i_sum = sums.begin();
      t_Column :: iterator i_2sum;
      for(; i_indiv != i_end; ++i_indiv, ++i_sum )
      {
        const t_Fitness& fit1  = Modifier::const_innermost( *i_indiv ).fitness();
        i_2indiv = i_indiv + 1;  
        i_2sum = i_sum + 1;
        --(*i_sum);
        // "Minimizes" the objectives/fitnesses by default
        for(; i_2indiv != i_end; ++i_2indiv, ++i_2sum )
          switch( fit1.compare( Modifier::const_innermost( *i_2indiv ).fitness() ) )
          {
            case Fitness::INDIFFERENT: break;
            case Fitness::WEAKER: --(*i_sum); break;
            case Fitness::STRONGER: --(*i_2sum); break;
            case Fitness::EQUAL: --(*i_sum); --(*i_2sum); break;
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
        const TiXmlElement *parent = &_node;

  // Don't need this if there is only one sharing function...
  // Also, it may be better to do everything as attributes directly
  //     if( name == "Sharing" )
  //       parent = &_node;
  //     else
  //       parent = _node.FirstChildElement("Sharing");
  //     
  //     for(; parent; parent = parent->NextSiblingElement("Sharing") )
  //     {
  //       if ( not parent->Attribute("type") ) continue;
  //       name = parent->Attribute("type");
  //       if(    name == "Triangular"
  //           or name == "triangular" ) break;
  //     }

        if ( not parent ) return false;
        if ( parent->Attribute("alpha") )
        {
          int d = 0;
          parent->Attribute("alpha", &d );
          __DOASSERT( d < 1, 
                         "Invalid alpha = " << d 
                      << " in Sharing::Triangular.\n" )
            
          alpha = types::t_unsigned( std::abs(d) );
        }
        if ( parent->Attribute("d0") )
        {
          double d = 0;
          parent->Attribute("d0", &d );
          __DOASSERT( Fuzzy::le(d, 0e0),
                         "Invalid d0 = " << d 
                      << " in Sharing::Triangular.\n" )
            
          d_0 = t_ScalarFitnessQuantity( d );
        }
        return distance.Load( _node );
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
         GeneralHamming<T_GAOPTRAITS>::operator()( const t_Individual &_i1, 
                                                   const t_Individual &_i2) const
         {
           __ASSERT( _i1.Object().Container().size() != _i1.Object().Container().size(),
                        "Cannot compute distance between individuals of different sizes\n"
                     << "Size of individual A: " << _i1.Object().Container().size() << "\n"
                     << "Size of individual B: " << _i2.Object().Container().size() << "\n")
           typename t_Object::const_iterator i_bit2 = _i2.Object().begin();
           typename t_Object::const_iterator i_bit1 = _i1.Object().begin();
           typename t_Object::const_iterator i_bit1_end = _i1.Object().end();
           t_ScalarFitnessQuantity result(0);
           for(; i_bit1 != i_bit1_end; ++i_bit1, ++i_bit2 )
             result += std::abs(   t_ScalarFitnessQuantity(*i_bit1) 
                                 - t_ScalarFitnessQuantity(*i_bit2) );
           return result;
         }

      template<class T_GAOPTRAITS>
      typename Hamming<T_GAOPTRAITS>::t_ScalarFitnessQuantity 
         Hamming<T_GAOPTRAITS>::operator()( const t_Individual &_i1, const t_Individual &_i2) const
         {
           __ASSERT( _i1.Object().Container().size() != _i1.Object().Container().size(),
                        "Cannot compute distance between individuals of different sizes\n"
                     << "Size of individual A: " << _i1.Object().Container().size() << "\n"
                     << "Size of individual B: " << _i2.Object().Container().size() << "\n")
           typedef typename t_GATraits::t_QuantityTraits::t_ScalarQuantityTraits t_SQTraits;
           typename t_Object::const_iterator i_bit2 = _i2.Object().begin();
           typename t_Object::const_iterator i_bit1 = _i1.Object().begin();
           typename t_Object::const_iterator i_bit1_end = _i1.Object().end();
           t_ScalarFitnessQuantity result(0);
           for(; i_bit1 != i_bit1_end; ++i_bit1, ++i_bit2 )
             result +=  t_ScalarFitnessQuantity(t_SQTraits::eq( *i_bit1, *i_bit2) ? 1: 0);
           return result;
         }

    } // namespace Distance


    template<class T_GATRAITS> 
    Base<T_GATRAITS>* new_from_xml( const TiXmlElement &_node, 
                                    typename T_GATRAITS :: t_Evaluator *_eval )
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
          container->push_back( new_Niche_from_xml<T_GATRAITS>( *parent, _eval ) );
        if(     T_GATRAITS::t_QuantityTraits::is_vector
            and ( name == "Pareto" or name == "pareto" ) )
          container->push_back( new ParetoRanking<T_GATRAITS>() );
        else if ( _eval )
          container->push_back( (Base<T_GATRAITS>*) _eval->Load_Scaling( _node ) );
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
    Base<T_GATRAITS>* new_Niche_from_xml( const TiXmlElement &_node,
                                          typename T_GATRAITS :: t_Evaluator *_eval )
    {
      if( not _node.Attribute("distance") ) return NULL;
      std::string name = _node.Attribute("distance");
      Base<T_GATRAITS> *result;
      const boost::regex gere("\\(G|g)eneral\\s*(H|h)amming");
      boost::match_results<std::string::const_iterator> what;
      if( boost::regex_search( name, what, gere ) )
      {
        typedef Niching< Sharing::Triangular<
                             Distance::GeneralHamming<T_GATRAITS> > > t_Niche;
        result = new t_Niche;
      }
      else if( name == "Hamming" or name == "hamming" )
      {
        typedef Niching< Sharing::Triangular< Distance::Hamming<T_GATRAITS> > > t_Niche;
        result = new t_Niche;
      }
      else if ( _eval ) result = ( Base<T_GATRAITS>* ) _eval->Load_Niche( _node );
       
      try { if ( result and result->Load(_node ) ) return result; }
      catch( std::exception &_e )
      {
        delete result; 
        __THROW_ERROR("Error while reading Niche input.\n"; )
      }
      if ( result ) delete result;
      return NULL;
    }

  } // namespace Scaling
} // namespace LaDa
#endif // _RANKING_IMPL_H_

//
//  Version: $Id$
//
#ifndef _RANKING_H_
#define _RANKING_H_

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
    class Base : public eoUF< typename T_GATRAITS::t_Population&, void> {};

  template<class T_GATRAITS> 
  Base<T_GATRAITS>* new_from_xml( const TiXmlElement &_node );

  template<class T_GATRAITS> 
  Base<T_GATRAITS>* new_Niche_from_xml( const TiXmlElement &_node );

  template<class T_GATRAITS>
  class Container : public Base<T_GATRAITS>
  {
    public:
      typedef T_GATRAITS t_GATraits;
    protected:
      typedef typename t_GATraits :: t_Population t_Population;
      typedef Base<t_GATraits> t_Base;
      typedef std::list<t_Base*> t_Container;

    protected:
      t_Container container;

    public:
      Container() {}
      Container( const Container &_c ) : container(_c.container) {}
      ~Container();

      void push_back( t_Base* _ptr) { if( _ptr ) container.push_back( _ptr ); }
      void push_front( t_Base* _ptr) { if( _ptr ) container.push_front( _ptr ); }

      void operator()( t_Population & _pop );

      types::t_unsigned size() const { return container.size(); }
  };


  template<class T_SHARING>
  class Niching : Base<typename T_SHARING::t_GATraits >
  {
    public:
      typedef T_SHARING t_Sharing;
      typedef typename t_Sharing :: t_GATraits t_GATraits;
    protected:
      typedef Base< t_GATraits > t_Base;
      typedef typename t_GATraits :: t_Individual t_Individual;
      typedef typename t_GATraits :: t_Fitness t_Fitness;
      typedef typename t_GATraits :: t_Population t_Population;
      typedef typename t_Fitness  :: t_Quantity t_FitnessQuantity;
      typedef typename t_GATraits :: t_ScalarFitness t_ScalarFitness;
      typedef typename t_ScalarFitness :: t_Quantity t_ScalarFitnessQuantity;
      typedef std::vector< t_ScalarFitnessQuantity >  t_Column;

    protected:
      t_Sharing sharing;
      t_Column sums;

    public:
  
      Niching() : t_Base(), sharing() {}
      Niching( const TiXmlElement &_node ) : t_Base(), sharing( _node ) {}
      Niching( const Niching &_n) : t_Base(_n), sharing(_n.sharing), sums(_n.sums) {}
      ~Niching() {}
  
      void operator()(t_Population& _pop);
      bool Load( const TiXmlElement &_node )
        { return sharing.Load(_node); }

    protected:
      void map(const t_Population &_pop);
  };

  template<class T_GATRAITS>
  class ParetoRanking : Base< T_GATRAITS >
  {
    public:
      typedef T_GATRAITS t_GATraits;
    protected:
      typedef Base< t_GATraits > t_Base;
      typedef typename t_GATraits :: t_Individual    t_Individual;
      typedef typename t_GATraits :: t_Population    t_Population;
      typedef typename t_GATraits :: t_Fitness       t_Fitness;
      typedef typename t_Fitness :: t_Quantity       t_FitnessQuantity;
      typedef typename t_Fitness :: t_ScalarFitness  t_ScalarFitness;
      typedef typename t_ScalarFitness :: t_Quantity t_ScalarFitnessQuantity;
      typedef std::vector< types::t_unsigned >       t_Column;

    protected:
      t_Column sums;

    public:
  
      ParetoRanking() : t_Base() {}
      ParetoRanking( const TiXmlElement &_node ) : t_Base() {}
      ParetoRanking( const ParetoRanking &_n) : t_Base(_n), sums(_n.sums) {}
      ~ParetoRanking() {}
  
      void operator()(t_Population& _pop);

      bool Load( const TiXmlElement &_node )  { return true; }
     
    protected:
      void map(const t_Population &_pop);
  };
} // namespace Ranking

namespace Distance
{
  template<class T_GATRAITS>
  class Hamming
  {
    public:
      typedef T_GATRAITS t_GATraits;
    protected:
      typedef typename t_GATraits :: t_Individual t_Individual;
      typedef typename t_GATraits :: t_Object t_Object;
      typedef typename t_GATraits :: t_ScalarFitness t_ScalarFitness;
      typedef typename t_ScalarFitness :: t_Quantity t_ScalarFitnessQuantity;

      t_ScalarFitnessQuantity operator()( const t_Individual &_i1, const t_Individual &_i2) const;
  };
} // namespace Distance

namespace Sharing
{
  template <class T_DISTANCE>
  class Triangular 
  {
    public:
      typedef T_DISTANCE t_Distance;
      typedef typename t_Distance :: t_GATraits t_GATraits;
    protected:
      typedef typename t_GATraits :: t_Individual t_Individual;
      typedef typename t_GATraits :: t_Fitness t_Fitness;
      typedef typename t_GATraits :: t_ScalarFitness t_ScalarFitness;
      typedef typename t_ScalarFitness :: t_Quantity t_ScalarFitnessQuantity;

    protected:
      t_Distance distance;
      types::t_unsigned alpha;
      t_ScalarFitnessQuantity sigma;

    public:
      Triangular() : distance(), alpha(1), sigma(1)  {}
      Triangular  ( const TiXmlElement &_node ) 
                 : distance(), alpha(1), sigma(1) { Load(_node); }
      Triangular   (const Triangular &_t )
                 : distance(_t.distance), alpha(_t.alpha), 
                   sigma(_t.sigma)  {}
      ~Triangular() {}

      bool Load( const TiXmlElement &_node );
      t_ScalarFitnessQuantity operator()(const t_Individual& _i1, const t_Individual& _i2 ) const;

      std::string print_out() const;
      bool is_valid() const { return distance != NULL; }
  };  
}

#include "ranking.impl.h"


#endif // _RANKING_H_

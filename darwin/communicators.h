//
//  Version: $Id$
//
#ifndef  _DARWIN_COMMUNICATORS_H_
#define  _DARWIN_COMMUNICATORS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <opt/types.h>
#include <opt/debug.h>
#include <mpi/mpi_object.h>

namespace GA
{
  namespace mpi 
  {

    class Graph : ::mpi::Base
    {
      public:
        template<class T_GATRAITS> class Breeder;
      protected: 
        MPI::IntraComm *herd_comm;
        MPI::IntraComm *farmer_comm;
        MPI::GraphComm  *graph_comm;
        types::t_unsigned pools;
        //! Type of the process
        enum t_Type 
        {
          COLLECTIVIST, //!<  Equivalent to pool == 1.
          FARMER, //!< The head boss.
          FARMHANDS, //!< Will do nothing.
          BULL,   //!< Takes orders from the farmer.
          COW     //!< Takes orders from one specific bull
        }
        t_Type type;

      public:  
        //! \brief Constructor and Initializer.
        //! \details The MPI::WORLD_COMM is duplicated and the duplicated is
        //!          pointed to by ::mpi::Base::comm.
        Graph() : ::mpi::Base::Base(), herd_comm(NULL), farmer_comm(NULL),
                  graph_comm(NULL), pools(0), type(FARMHANDS)
           { comm = &_comm->Clone(); }
        Graph   ( MPI::Comm &_comm )
              : ::mpi::Base::Base(_comm), herd_comm(NULL), farmer_comm(NULL),
                graph_comm(NULL), pools(0), type(FARMHANDS)
          { comm = &_comm->Clone(); }
        ~Graph();

        template< class T_CONDITION > void init( T_CONDITION &_condition );

        template< class T_GATRAITS >
        Breeder<T_GATRAITS>* create_breeder( eoState &_state );

      protected:
        template< class T_CONDITION > void init_collectivists( T_CONDITION &_condition );
    };

    template<class T_GATRAITS>
    class Graph::Breeder : public eoBreed<typename T_GATRAITS::t_Individual>
    {
      public:
        typedef T_GATRAITS t_GATraits; //!< all %GA traits
      protected:
        //! type of an individual
        typedef typename t_GATraits::t_Individual  t_Individual; 

      protected:
        //! A selection operator for the obtention of parents
        eoSelectOne<t_Individual>* select;
        //! Mating operator
        eoGenOp<t_Individual> *op;
        //! Generation counter
        GenCount *age;
        //! Number of offspring to change
        eoHowMany howMany;
        //! mpi topology
        Graph *topo;

      public:
        //! Constructor and Initializer
        Breeder   ( Graph *_topo )
                : select(NULL), op(NULL), age(NULL),
                  howmany(0), topo(_topo) {}

        //! Copy Constructor
        Breeder ( Breeder<t_Individual> & _breeder )
                : select( _breeder.topo ), op(_breeder.op), age(_breeder.age),
                  howmany(_breeder.howMany), topo( _breeder.topo ) {}

        //! Sets the selector.
        void set( eoSelectOne<t_Individual> *_select ){ select = _select; }
        //! Sets the breeding operators
        void set( eoSelectOne<t_Individual> *_op ){ op = _op; }
        //! Sets the replacement rate
        void set( types::t_real _rep ){ howMany = _rep; }
        //! Sets the generation counter.
        void set( GenCount *_age ){ age = _age; }
        //! Sets the topology.
        void set( Graph *_topo) { topo(_topo); }
        //! Destructor
        virtual ~Breeder() {}
        //! EO required.
        virtual std::string className() const { return "GA::mpi::Graph::Breeder"; }
    }

    template<class T_GATRAITS>
    class FarmhandGraphBreeder : public Graph::Breeder<T_GATRAITS>
    {
      public:
        typedef T_GATRAITS t_GATraits; //!< all %GA traits

      protected:
        //! all individual traits
        typedef typename t_GATraits::t_IndivTraits t_IndivTraits;
        //! type of an individual
        typedef typename t_GATraits::t_Individual  t_Individual; 
        //! type of the population
        typedef typename t_GATraits::t_Population  t_Population; 
        //! Base class type.
        typedef Graph::Breeder<T_GATRAITS> t_Base;

      public:
        //! Constructor.
        FarmhandGraphBreeder( Graph *_topo ) : t_Base( _topo );

        //! Creates \a _offspring population from \a _parent
        void operator()(const t_Population& _parents, t_Population& _offspring)
          { return; }
   
        ///! The class name. EO required
        virtual std::string className() const { return "GA::mpi::FarmhandGraphBreeder"; }
    };

    template<class T_GATRAITS>
    class FarmerGraphBreeder : public Graph::Breeder<T_GATRAITS>
    {
      public:
        typedef T_GATRAITS t_GATraits; //!< all %GA traits

      protected:
        //! all individual traits
        typedef typename t_GATraits::t_IndivTraits t_IndivTraits;
        //! type of an individual
        typedef typename t_GATraits::t_Individual  t_Individual; 
        //! type of the population
        typedef typename t_GATraits::t_Population  t_Population; 
        //! Base class type.
        typedef Graph::Breeder<T_GATRAITS> t_Base;

      public:
        //! Constructor.
        FarmerGraphBreeder( Graph *_topo ) : t_Base( _topo );

        //! Creates \a _offspring population from \a _parent
        void operator()(const t_Population& _parents, t_Population& _offspring);
   
        ///! The class name. EO required
        virtual std::string className() const { return "GA::mpi::FarmerGraphBreeder"; }
    };
  } // namespace mpi
} // namespace GA

#include "comminucators.impl.h"

#endif

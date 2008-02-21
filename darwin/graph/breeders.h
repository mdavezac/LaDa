//
//  Version: $Id$
//
#ifndef  _GRAPH_BRREEDERS_H_
#define  _GRAPH_BRREEDERS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef _MPI

#include <opt/types.h>
#include <opt/debug.h>
#include <mpi/mpi_object.h>
#include <darwin/breeder.h>
#include "comm.h"

namespace GA
{
  namespace mpi 
  {


    namespace Graph
    {
      //! Contains all Breeder related stuff in the mpi::Graph Topology.
      namespace Breeder
      {

        //! \brief Base class for breeders in the GA::mpi::Graph topology
        //! \details Contains helper functions which may be of help to any of the
        //!          specific breeders.
        template<class T_GATRAITS>
        class Base : public eoBreed<typename T_GATRAITS::t_Individual>
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
    
          public:
            //! Constructor and Initializer
            Base () : select(NULL), op(NULL), age(NULL),
                      howmany(0), topo(_topo) {}
            //! Copy Constructor
            Base   ( Base<t_Individual> & _breeder )
                 : select( _breeder.topo ), op(_breeder.op), age(_breeder.age),
                   howmany(_breeder.howMany), topo( _breeder.topo ) {}
            //! Destructor
            virtual ~Base() {}
    
            //! Sets the selector.
            void set( eoSelectOne<t_Individual> *_select ){ select = _select; }
            //! Sets the breeding operators
            void set( eoSelectOne<t_Individual> *_op ){ op = _op; }
            //! Sets the replacement rate
            void set( types::t_real _rep ){ howMany = _rep; }
            //! Sets the generation counter.
            void set( GenCount *_age ){ age = _age; }
            //! Sets the topology.
            void set( Topology *_topo) { topo(_topo); }
            //! EO required.
            virtual std::string className() const
              { return "GA::mpi::Graph::Breeder::Base"; }
        }
    
        //! A breeder class which does nothing.
        template<class T_GATRAITS>
        class Farmhand : public Base<T_GATRAITS>
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
            typedef Base<t_GATraits> t_Base;
    
          public:
            //! Constructor.
            Farmhand( Topology *_topo ) : t_Base( _topo );
    
            //! Creates \a _offspring population from \a _parent
            void operator()(const t_Population& _parents, t_Population& _offspring) {};
       
            ///! The class name. EO required
            virtual std::string className() const
              { return "GA::mpi::Graph::Breeder::Farmhand"; }
        };
    
        //! \brief A breeder class to rule them all.
        //! \details This functor dispatches commands to the bulls, such as breed
        //!          one and stop breeding. 
        template<class T_GATRAITS>
        class Farmer : private Comm::Farmer< Farmer<T_GATRAITS> >, 
                       public Base<T_GATRAITS>
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
            typedef Base<t_GATraits> t_Base;
            //! Communication base class
            typedef Comm::Farmer< Farmer<t_GATraits> > t_CommBase;
    
          protected:
            types::t_unsigned target;
            t_Population *offspring;
            
    
          public:
            //! Constructor.
            Farmer   ( Topology *_topo )
                   : t_CommBase( _topo ), t_Base(),
                     target(0), offspring(NULL) {};
    
            //! Creates \a _offspring population from \a _parent
            void operator()(const t_Population& _parents, t_Population& _offspring);
       
            //! The class name. EO required
            virtual std::string className() const
              { return "GA::mpi::Graph::Breeder::Farmer"; }
    
            //! \brief Response to WAITING request.
            //! \details Receives the individual. If still more are to be
            //!          created, sends t_CommBase::t_Commands::DONE, otherwise
            //!          t_CommBase::t_Commands::GO. The persistent request is
            //!          only activated in the latter case.
            void onWait( types::t_int _bull );
        };
    
        template<class T_GATRAITS>
        class Bull : private Comm::Bull< Bull<T_GATRAITS> >, public Base<T_GATRAITS>
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
            typedef Base<t_GATraits> t_Base;
            //! Communication base class
            typedef Comm::Bull< Bull<t_GATraits> > t_CommBase;
    
            //! Tag for communications with the cows
            const MPI::INT COWTAG = 2;
    
    
          public:
            //! Constructor.
            Bull   ( Topology *_topo )
                 : t_CommBase(_topo ), t_Base() {}
    
            //! Creates \a _offspring population from \a _parent
            void operator()(const t_Population& _parents, t_Population& _offspring);
       
            //! The class name. EO required
            virtual std::string className() const
              { return "GA::mpi::Graph::Breeder::Bull"; }
        };
    
        //! Just a typedef for Comm::BaseCow.
        template< class T_GATRAITS >
        class Cow : public Comm::LaNormande< T_GATRAITS, Base > {};
    
  } // namespace mpi
} // namespace GA

#include "breeders.impl.h"

#endif
#endif

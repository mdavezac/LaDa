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


    namespace Graph
    {
      //! \brief Contains all evalaution related stuff in the mpi::Graph Topology.
      //! \details The classes below are meant to be used as derived from
      //!          "standard" Evaluation classes. They will simply overwrite
      //!          the behavior of the Evaluation::Base::operator()(
      //!          t_Population &_parent, t_Population &_offspring ) routine.
      //!          Instead of the standard approach to evaluating offsprings,
      //!          each invalid individual is dispatched by the farmer to the
      //!          herds for evaluation. Individuals are dispatched as herds
      //!          become available. The classes from which the following
      //!          classes are derive is defined through a template. These
      //!          classes also derive from the Comm classes, \e via the CRT
      //!          mechanism.
      namespace Evaluation
      {

        //! An evaluation class which does nothing.
        template<class T_BASE>
        class Farmhand
        {
          public:
            //! Base class type
            typedef T_BASE t_Base;
    
          protected:
            //! all %GA traits
            typedef typename t_Base::t_GATraits t_GATraits;
            //! all individual traits
            typedef typename t_GATraits::t_IndivTraits t_IndivTraits;
            //! type of the population
            typedef typename t_GATraits::t_Population  t_Population; 
    
          public:
            //! Constructor.
            Farmhand();
    
            //! Creates \a _offspring population from \a _parent
            void operator()(const t_Population& _parents, t_Population& _offspring) {};
        };
    
        //! \brief A breeder class to rule them all.
        //! \details This functor dispatches commands to the bulls. More
        //!          specifically, it makes a list of unknown individuals and
        //!          dispatches them for evaluation to the bulls.
        template<class T_BASE>
        class Farmer : private Comm::Farmer< Farmer >, public T_BASE
        {
          public:
            //! Base class type
            typedef T_BASE t_Base;
    
          protected:
            //! all %GA traits
            typedef typename t_Base::t_GATraits t_GATraits;
            //! all individual traits
            typedef typename t_GATraits::t_IndivTraits t_IndivTraits;
            //! type of an individual
            typedef typename t_GATraits::t_Individual  t_Individual; 
            //! type of the population
            typedef typename t_GATraits::t_Population  t_Population; 
            //! Communication base class
            typedef Comm::Farmer< t_GATraits, Farmer > t_CommBase;
    
          protected:
            types::t_unsigned target;
            t_Population *offspring;
            std::list< t_Population :: iterator > unknowns;
    
          public:
            //! Constructor.
            Farmer   ( Topology *_topo )
                   : t_CommBase(_topo->Com() ), t_Base( _topo->Comm() ),
                     target(0), offspring(NULL) {};
    
            //! Creates \a _offspring population from \a _parent
            void operator()(const t_Population& _parents, t_Population& _offspring);
       
            //! The class name. EO required
            virtual std::string className() const { return "GA::mpi::Graph::Breeder::Farmer"; }
    
            // Response to WAITING request
            void onWait( types::t_int _bull );
            // Response to REQUESTINGTABOOCHECK request
            void onTaboo( types::t_int _bull );
            // Response to REQUESTINGOBJECTIVE request
            void onObjective( types::t_int _bull );
            // Response to REQUESTINGHISTORYCHECK request
            void onHistory( types::t_int _bull );
        };
    
        template<class T_GATRAITS>
        class Bull : private Comm::Bull<T_GATRAITS, Bull>, public Base<T_GATRAITS>
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
            typedef Base<T_GATRAITS> t_Base;
            //! Communication base class
            typedef Comm::Bull< T_GATRAITS, Farmer> t_CommBase;
    
            //! Tag for communications with the cows
            const MPI::INT COWTAG = 2;
    
    
          public:
            //! Constructor.
            Bull   ( Topology *_topo )
                 : t_CommBase(_topo->Com() ), t_Base( _topo->Comm() ) {}
    
            //! Creates \a _offspring population from \a _parent
            void operator()(const t_Population& _parents, t_Population& _offspring);
       
            //! The class name. EO required
            virtual std::string className() const { return "GA::mpi::Graph::Breeder::Bull"; }
        };
    
        template<class T_GATRAITS>
        class Cow : private CommBull< T_GATRAITS, Cow >, public Base<T_GATRAITS>
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
            //! Communication base class
            typedef Comm::Cow< T_GATRAITS, Cow > t_CommBase;
    
            const MPI::INT TAG = 2;
    
    
          public:
            //! Constructor.
            Cow   ( Topology *_topo )
                 : t_CommBase(_topo->Com() ), t_Base( _topo->Comm() ) {}
    
            //! Creates \a _offspring population from \a _parent
            void operator()(const t_Population& _parents, t_Population& _offspring);
              { while( t_CommBase :: wait() != t_CommBase::DONE ); }
       
            //! The class name. EO required
            virtual std::string className() const { return "GA::mpi::Graph::Breeder::Cow"; }
        };
      } // namespace Evaluation
    } // namespace Graph
  } // namespace mpi
} // namespace GA

#include "comminucators.impl.h"

#endif

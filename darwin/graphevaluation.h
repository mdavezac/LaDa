//
//  Version: $Id$
//
#ifndef  _DARWIN_COMMUNICATORS_H_
#define  _DARWIN_COMMUNICATORS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <list>
#include <pair>

#include <opt/types.h>
#include <opt/debug.h>
#include <mpi/mpi_object.h>
#include "evaluation.h"

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
        class Farmer : private Comm::Farmer< Farmer<T_BASE> >, public T_BASE
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
            //! Type of the history
            typedef GA::History<t_Individual> t_History;
            //! Communication base class
            typedef Comm::Farmer< Farmer<T_BASE> > t_CommBase;
            //! first holds iterator, second holds assigned proc.
            typedef std::pair< typename t_Population :: iterator,
                               types::t_int > t_Unknown;
            //! Container of individuals needing evaluation.
            typedef std::list< t_Unknown > t_Unknowns;
    
          protected:
            types::t_unsigned target;
            t_Population *offspring;
            t_Unknowns unknowns;
    
          public:
            //! Constructor.
            Farmer   ( Topology *_topo )
                   : t_CommBase( _topo ), T_BASE()
                     target(0), offspring(NULL) {};
    
            //! Creates \a _offspring population from \a _parent
            void operator()(const t_Population& _parents, t_Population& _offspring);
    
            //! The class name. EO required
            virtual std::string className() const
              { return "GA::mpi::Graph::Evaluation::Farmer"; }

            //! Sets taboo pointer
            set( Taboo_Base<t_Individual> *_taboo )
              { t_CommBase :: taboos = t_Base :: taboos = _taboos; }
            //! Sets objective pointer
            set(  typename t_ObjectiveType::Vector*  _o )
              { t_CommBase :: objective = t_Base :: objective = _o; }
            //! Sets objective pointer
            set(  typename t_Store::Base*  _s )
              { t_CommBase :: store = t_Base :: store =  _s; }
            //! Sets history pointer
            set(  typename t_History*  _h)
              { t_CommBase :: history = _h; set_Base_history<t_Base>(_h); }

          protected:
            //! Response to WAITING request
            void onWait( types::t_int _bull );
            //! Sets history for t_Base != Evaluation::WithHistory.
            template < class TT_BASE > set_base_history( t_History *_history ) {}
        };
    
        template<class T_BASE>
        class Bull : private Comm::Bull< Bull<T_BASE> >, public T_BASE
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
            virtual std::string className() const
              { return "GA::mpi::Graph::Evaluation::Bull"; }
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
            virtual std::string className() const
              { return "GA::mpi::Graph::Breeder::Cow"; }
        };
      } // namespace Evaluation
    } // namespace Graph
  } // namespace mpi
} // namespace GA

#include "comminucators.impl.h"

#endif

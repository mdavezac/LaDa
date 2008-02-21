//
//  Version: $Id$
//
#ifndef  _GRAPH_EVALUATION_H_
#define  _GRAPH_EVALUATION_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef _MPI

#include <list>
#include <pair>

#include <opt/types.h>
#include <opt/debug.h>
#include <mpi/mpi_object.h>
#include "evaluation.h"
#include "graphevaluation.h"

namespace GA
{
  namespace mpi 
  {


    namespace Graph
    {
      //! \brief Contains all evaluation related stuff in the mpi::Graph Topology.
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
        template< class T_GATRAITS, template < class > class T_BASE>
        class Farmer : private Comm::Farmer< Farmer<T_GATRAITS, T_BASE> >,
                       public T_BASE< GA::Traits< Evaluator::Farmer<T_GATRAITS>,
                                                  typename T_GATRAITS :: t_Population,
                                                  typename T_GATRAITS :: t_Islands > >
        {
          typedef Traits::GA< Evaluator<T_GATRAITS>,
                              typename T_GATRAITS :: t_Population,
                              typename T_GATRAITS :: t_Islands > t_BaseTraits;
          public:
            //! Base class type with cache-evaluator.
            typedef T_BASE<t_BaseTraits> t_Base;
            //! all %GA traits
            typedef typename T_GATRAITS t_GATraits;
    
          protected:
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
            //! Type of the meta-evaluator.
            typedef GA::mpi::Graph::Evaluator::Farmer< t_GATraits > t_CacheEvaluator;
    
          protected:
            t_CacheEvaluator cache_eval;
    
          public:
            //! Constructor.
            Farmer   ( Topology *_topo )
                   : t_CommBase( _topo ), t_Base()
                     target(0), offspring(NULL) { t_Base::evaluator = &cache_eval; };
    
            //! Creates \a _offspring population from \a _parent
            void operator()(const t_Population& _parents, t_Population& _offspring);
    
            //! The class name. EO required
            virtual std::string className() const
              { return "GA::mpi::Graph::Evaluation::Farmer"; }

            //! Does nothing.
            set( t_Evaluator *_e ) {}
            //! Sets taboo pointer
            set( Taboo_Base<t_Individual> *_taboo )
              { t_CommBase :: taboos = t_Base :: taboos = _taboos; }
            //! Sets objective pointer
            set(  typename t_ObjectiveType::Vector*  _o )
              { t_CommBase :: objective = t_Base :: objective = _o; }
            //! Sets objective pointer
            set(  typename t_Store::Base*  _s )
              { t_CommBase :: store = t_Base :: store =  _s; }

          protected:
            //! Response to WAITING request
            void onWait( types::t_int _bull );
        };
    
        template< class T_GATRAITS, template < class > class T_BASE>
        class Bull : private Comm::Bull< Bull<T_GATRAITS, T_BASE> >,
                     public T_BASE< GA::Traits< Evaluator::Bull<T_GATRAITS>,
                                                typename T_GATRAITS :: t_Population,
                                                typename T_GATRAITS :: t_Islands > >
        {
          typedef Traits::GA< Evaluator<T_GATRAITS>,
                              typename T_GATRAITS :: t_Population,
                              typename T_GATRAITS :: t_Islands > t_BaseTraits;
          public:
            //! Base class type with meta-evaluator
            typedef T_BASE<t_BaseTraits> t_Base;
            //! all %GA traits
            typedef typename T_GATRAITS t_GATraits;
    
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
                 : t_CommBase( _topo ) {}
    
            //! Creates \a _offspring population from \a _parent
            void operator()(const t_Population& _parents, t_Population& _offspring);
       
            //! The class name. EO required
            virtual std::string className() const
              { return "GA::mpi::Graph::Evaluation::Bull"; }
        };

        
        //! Just a typedef for Comm::BaseCow.
        template< class T_GATRAITS, template < class > class T_BASE>
        class Cow : public Comm::LaNormande< T_GATRAITS, T_BASE > {};
    
      } // namespace Evaluation
    } // namespace Graph
  } // namespace mpi
} // namespace GA

#include "evaluation.impl.h"

#endif
#endif

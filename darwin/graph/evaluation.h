//
//  Version: $Id$
//
#ifndef  _GRAPH_EVALUATION_H_
#define  _GRAPH_EVALUATION_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <darwin/evaluation.h>

#ifdef _MPI

#include <opt/types.h>
#include <opt/debug.h>
#include <mpi/mpi_object.h>
#include <darwin/gatraits.h>

#include "evaluator.h"

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
            Farmhand() {}
    
            //! Creates \a _offspring population from \a _parent
            void operator()(const t_Population& _parents, t_Population& _offspring) {};
        };
    
        //! \brief A breeder class to rule them all.
        //! \details This functor dispatches commands to the bulls. More
        //!          specifically, it makes a list of unknown individuals and
        //!          dispatches them for evaluation to the bulls.
        template< class T_GATRAITS, template < class > class T_BASE>
        class Farmer : protected Comm::Farmer< T_GATRAITS,
                                               Farmer<T_GATRAITS, T_BASE> >,
                       public T_BASE< Traits::GA< Evaluator::Farmer<T_GATRAITS>,
                                                  typename T_GATRAITS :: t_Population,
                                                  typename T_GATRAITS :: t_Islands > >
        {
          friend class Comm::Farmer< T_GATRAITS, Farmer<T_GATRAITS, T_BASE> >;
          typedef Traits::GA< Evaluator::Farmer<T_GATRAITS>,
                              typename T_GATRAITS :: t_Population,
                              typename T_GATRAITS :: t_Islands > t_BaseTraits;
          public:
            //! Base class type with cache-evaluator.
            typedef T_BASE<t_BaseTraits> t_Base;
            //! all %GA traits
            typedef T_GATRAITS t_GATraits;
    
          protected:
            //! Type of this class.
            typedef Farmer< t_GATraits, T_BASE >         t_This;
            //! Type of the original evaluator.
            typedef typename t_GATraits::t_Evaluator     t_Evaluator;
            //! all individual traits.
            typedef typename t_GATraits::t_IndivTraits   t_IndivTraits;
            //! type of an individual.
            typedef typename t_GATraits::t_Individual    t_Individual; 
            //! type of the population.
            typedef typename t_GATraits::t_Population    t_Population; 
            //! Type of the history.
            typedef GA::History<t_Individual>            t_History;
            //! Communication base class.
            typedef Comm::Farmer< T_GATRAITS, t_This >   t_CommBase;
            //! Type of the meta-evaluator.
            typedef typename t_BaseTraits :: t_Evaluator t_CacheEvaluator;
    
          protected:
            //! Meta-evaluator for distributing evaluations to procs.
            t_CacheEvaluator cache_eval;
    
          public:
            //! Constructor.
            Farmer   ( Topology *_topo )
                   : t_CommBase( _topo ), t_Base(), cache_eval() 
                     { t_Base::evaluator = &cache_eval; };
    
            //! Creates \a _offspring population from \a _parent
            void operator()(t_Population& _parents, t_Population& _offspring);
    
            //! The class name. EO required
            virtual std::string className() const
              { return "GA::mpi::Graph::Evaluation::Farmer"; }

            //! Does nothing.
            void set( t_Evaluator *_e ) {}
            //! Sets taboo pointer
            void set( GA::Taboo_Base<t_Individual> *_taboo )
              { t_CommBase :: taboos = t_Base :: taboos = _taboo; }
            //! Sets objective pointer
            void set(  typename GA::Objective::Types<t_GATraits>::Vector*  _o )
              { t_CommBase :: objective = t_Base :: objective = _o; }
            //! Sets objective pointer
            void set(  typename GA::Store::Base<t_GATraits>*  _s )
              { t_CommBase :: store = t_Base :: store =  _s; }


          protected:
            //! Response to WAITING request
            void onWait( types::t_unsigned _bull );
            using t_Base::evaluate;
        };
    
        template< class T_GATRAITS, template < class > class T_BASE>
        class Bull : protected Comm::Bull< T_GATRAITS, Bull<T_GATRAITS, T_BASE> >,
                     public T_BASE< Traits::GA< Evaluator::Bull<T_GATRAITS>,
                                                typename T_GATRAITS :: t_Population,
                                                typename T_GATRAITS :: t_Islands > >
        {
          typedef Traits::GA< Evaluator::Bull<T_GATRAITS>,
                              typename T_GATRAITS :: t_Population,
                              typename T_GATRAITS :: t_Islands > t_BaseTraits;
          public:
            //! Base class type with meta-evaluator
            typedef T_BASE<t_BaseTraits> t_Base;
            //! all %GA traits
            typedef T_GATRAITS t_GATraits;
    
          private:
            //! This type.
            typedef Bull<t_GATraits, T_BASE >            t_This;
            //! all individual traits.
            typedef typename t_GATraits::t_IndivTraits   t_IndivTraits;
            //! type of an individual.
            typedef typename t_GATraits::t_Individual    t_Individual; 
            //! type of the population.
            typedef typename t_GATraits::t_Population    t_Population; 
            //! Communication base class.
            typedef Comm::Bull< t_GATraits, t_This >     t_CommBase;
            //! Type of the original evaluator.
            typedef typename t_GATraits::t_Evaluator     t_Evaluator;
            //! Type of the meta-evaluator.
            typedef typename t_BaseTraits :: t_Evaluator t_BullEvaluator;
    
          protected:
            t_BullEvaluator metaeval;

          public:
            //! Constructor.
            Bull   ( Topology *_topo )
                 : t_CommBase( _topo ), t_Base(),
                   metaeval(_topo) { t_Base::evaluator = &metaeval; }
    
            //! Creates \a _offspring population from \a _parent
            virtual void operator()(t_Population& _parents, t_Population& _offspring);
       
            //! The class name. EO required
            virtual std::string className() const
              { return "GA::mpi::Graph::Evaluation::Bull"; }
            //! Does nothing.
            void set( t_Evaluator *_e ) { metaeval = *_e; }
        };

        
        //! Just a typedef for Comm::BaseCow.
        template< class T_GATRAITS, template < class > class T_BASE>
        class Cow : public Comm::LaNormande< T_GATRAITS, T_BASE > 
        {
          public:  
            //! All %GA types.
            typedef T_GATRAITS t_GATraits;
            //! Base class of LaNormande.
            typedef T_BASE<t_GATraits> t_Base;
            //! Constructor.
            Cow   ( Topology * _topo )
                : Comm::LaNormande< t_GATraits, T_BASE >( _topo ){}
        };
    
      } // namespace Evaluation
    } // namespace Graph
  } // namespace mpi
} // namespace GA

#include "evaluation.impl.h"

#endif
#endif

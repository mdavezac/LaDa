//
//  Version: $Id$
//
#ifndef  _GRAPH_EVALUATOR_H_
#define  _GRAPH_EVALUATOR_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef _MPI

#include <list>

#include <opt/types.h>
#include <opt/debug.h>
#include <mpi/mpi_object.h>
#include <darwin/evaluator.h>

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
      namespace Evaluator
      {

        //! \brief A meta-evaluator for bulls.
        //! \details It catches the evaluation calls and first broadcasts the
        //!          approapriate information to the bulls.
        template<class T_GATRAITS>
        class Bull : public GA::Evaluator< typename T_GATRAITS :: t_Individual >,
                     private Comm::Bull< T_GATRAITS, Bull<T_GATRAITS > >
        {
          public:
            //! all %GA traits
            typedef T_GATRAITS t_GATraits;
            //! type of an individual
            typedef typename t_GATraits::t_Individual           t_Individual; 
    
          protected:
            //! This class type 
            typedef Bull<t_GATraits>                            t_This;
            //! all individual traits.
            typedef typename t_GATraits::t_IndivTraits          t_IndivTraits;
            //! Base class type
            typedef GA::Evaluator<t_Individual>                 t_Base;
            //! Communication base class type
            typedef Comm::Bull<t_GATraits, t_This>              t_CommBase;
            //! type of the population
            typedef typename t_GATraits::t_Population           t_Population; 
            //! Type of the lamarckian traits, as declared in the base class
            typedef typename t_GATraits :: t_VA_Traits          t_VA_Traits;
            //! Type of the lamarckian gradients, as declared in the base class
            typedef typename t_VA_Traits :: t_QuantityGradients t_QuantityGradients;
            //! Type of the true evaluator
            typedef typename t_GATraits :: t_Evaluator          t_TrueEvaluator;

          protected:
            t_TrueEvaluator *trueEvaluator;
    
          public:
            //! Constructor.
            Bull ( Topology *_topo ) : t_Base(), t_CommBase( _topo ) {}
            //! Copy Constructor
            Bull ( const t_This &_c ) : t_Base( _c ) {}
    
            //! Commands herd to evaluate t_This::current_individual.
            void evaluate();
            //! Commands herd to perform a gradient evaluation.
            void evaluate_gradient( t_QuantityGradients& _grad );
            //! Commands herd to perform an evaluation with a gradient evaluation.
            void evaluate_with_gradient( t_QuantityGradients& _grad );
            //! Commands herd to perform the evaluation of one gradient.
            void evaluate_one_gradient( t_QuantityGradients& _grad,
                                        types::t_unsigned _pos);
       
            //! The class name. EO required
            virtual std::string className() const
              { return "GA::mpi::Graph::Evaluation::Bull"; }
            void set( t_TrueEvaluator *_eval )
            {
              if( not _eval ) Print :: out << "Evaluator pointer is null" << Print::endl;
              trueEvaluator = _eval; 
            }
            void init( t_Individual &_indiv )
            {
              __ASSERT( not trueEvaluator, "True Evaluator pointer not set" )
              t_Base::init( _indiv );
              trueEvaluator->init( _indiv );
            } 
        };
    
      } // namespace Evaluation
    } // namespace Graph
  } // namespace mpi
} // namespace GA

#include "evaluator.impl.h"

#endif
#endif
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
      namespace Evaluator

        //! \brief A meta-evaluator class for GA::mpi::Graph::Evaluation::Farmer.
        //! \details This class catches calls to GA::Evaluator::evaluate() and
        //           caches the individuals to be evaluated (by the bull
        //           processes). Farmer::onWait() erases the current individual
        //           from the unknown list, and returns the next individual to
        //           evaluate. It returns null if there are no more individuals
        //           to evaluate.
        template<class T_GATRAITS>
        class Farmer : public GA::Evaluator< typename T_GATRAITS :: t_Individual >
        {
          public:
            //! all %GA traits
            typedef typename T_GATRAITS t_GATraits;
    
          protected:
            //! This class type 
            typedef Farmer<t_GATraits> t_This;
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
            //! first holds pointer, second holds assigned proc.
            typedef std::pair< t_Individual*, types::t_int > t_Unknown;
            //! Container of individuals needing evaluation.
            typedef std::list< t_Unknown > t_Unknowns;
            //! Base class type
            typedef GA::Evaluator<t_Individual> t_Base;
    
          protected:
            t_Unknowns unknowns;
    
          public:
            //! Constructor.
            Farmer () : t_Base() {};
            //! Copy Constructor
            Farmer ( const t_This &_c ) : t_Base( _c ), unknowns( _c.unknowns ) {};

            //! Shoves individual which need evaluating in a list.
            void evaluate() 
             { unknowns.push_back( t_Unknown( &Modifier::innermost(i_indiv), -1 ) ); }

            //! Returns next unknown to evaluate.
            t_Individual* onWait( types::t_unsigned _bull );
            //! Should not be called: throws.
            void evaluate_gradient( t_QuantityGradients& _grad )
              { __THROW_ERROR( "Should not be called.\n" ) }
            //! Should not be called: throws.
            void evaluate_with_gradient( t_QuantityGradients& _grad )
              { __THROW_ERROR( "Should not be called.\n" ) }
            //! Should not be called: throws.
            void evaluate_one_gradient( t_QuantityGradients& _grad, types::t_unsigned _pos) 
              { __THROW_ERROR( "Should not be called.\n" ) }
        };
    
        template<class T_GATRAITS>
        class Bull : public GA::Evaluator< typename T_GATRAITS :: t_Evaluator >
        {
          public:
            //! all %GA traits
            typedef typename T_GATRAITS t_GATraits;
    
          protected:
            //! This class type 
            typedef Bull<t_GATraits> t_This;
            //! Base class type
            typedef GA::Evaluator<typename t_GATraits :: t_Evaluator> t_Base;
            //! type of an individual
            typedef typename t_Base::t_Individual  t_Individual; 
            //! all individual traits
            typedef typename t_Individual::t_IndivTraits t_IndivTraits;
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
            Bull () : t_Base() {}
            //! Copy Constructor
            Farmer ( const t_This &_c ) : t_Base( _c ) {}
    
            //! Commands herd to evaluate t_This::current_individual.
            void evaluate();
            //! Commands herd to perform a gradient evaluation.
            void evaluate_gradient( t_QuantityGradients& _grad );
            //! Commands herd to perform an evaluation with a gradient evaluation.
            void evaluate_with_gradient( t_QuantityGradients& _grad );
            //! Commands herd to perform the evaluation of one gradient.
            void evaluate_one_gradient( t_QuantityGradients& _grad, types::t_unsigned _pos);
       
            //! The class name. EO required
            virtual std::string className() const
              { return "GA::mpi::Graph::Evaluation::Bull"; }
        };
    
      } // namespace Evaluation
    } // namespace Graph
  } // namespace mpi
} // namespace GA

#include "graphevaluator.impl.h"

#endif

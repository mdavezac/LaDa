//
//  Version: $Id$
//
#ifndef _MULTIOB_EVALUATION_H_
#define _MULTIOB_EVALUATION_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <list> 

#include <eo/eoPop.h>

#include <tinyxml/tinyxml.h>

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif
#include "print/stdout.h"

#include "objective.h"
#include "store.h"
#include "taboos.h"
#include "gatraits.h"

/** \ingroup Genetic 
 * @{*/
//! \brief Contains wrapper classes which make call to an Evaluator's functional routines
//! \details The classes in this namespace form the funnel through which all
//! requests from %GA to evaluate any object are sent. The goal is to capture any
//! such request for storage puposes <em>via</em> the classes in the Store
//! namespace, and for taboo (history) purposes <em>via</em> Taboo::History. It
//! also allows us to separate the functional returned by the functional and
//! stored in Base::t_Individual::quantities() from the fitness of an individual,
//! implemented in Fitness and kept in Individual %type objects. Fitness and
//! quantities are linked <em>via</em> the classes implemented in Objective.
//! At this point, there are two classes:
//!  - Base: declares a simple base class with a complete set of storage behaviors
//!  - WithHistory: can also keep tab on all evaluated objects past and present
//! through a Taboo::History object. 
//!  .
//! \sa Store, Store::Base, Taboo::History
namespace Evaluation 
{ 
  //! \brief Funnel through which all calls from %GA to the %Evaluator's
  //! functional routines are routed.
  //! \details It implements obvious routines such as Base::evaluate(), and
  //! Base::evaluate_gradient, but also a less obvious though necessary storage
  //! ability <em>via</em> through a referenced to a Store::Base objet.
template< class T_GATRAITS >
  class Base
  {
    public:
      //! Should contain all GA definitions \sa Traits::GA
      typedef T_GATRAITS  t_GATraits;

    protected:
      typedef typename t_GATraits :: t_Individual                  t_Individual;
      typedef typename t_GATraits :: t_Evaluator                   t_Evaluator;
      typedef typename t_GATraits :: t_QuantityTraits              t_QuantityTraits;
      typedef typename t_QuantityTraits :: t_Quantity              t_Quantity;
      typedef typename t_QuantityTraits :: t_ScalarQuantity        t_ScalarQuantity;
      typedef typename t_GATraits :: t_Population                  t_Population;
      typedef typename Objective :: Types < t_GATraits > :: Vector t_Objective;
      typedef Store :: Base<t_GATraits>                            t_Store;
      typedef typename t_GATraits :: t_VA_Traits                   t_VATraits;
      typedef typename t_VATraits :: t_QuantityGradients           t_QuantityGradients;
      typedef typename t_VATraits :: t_Type                        t_VA_Type;
      typedef typename t_GATraits :: t_Fitness :: t_ScalarFitness  t_Fitness;
      typedef typename t_Fitness :: t_Quantity                     t_FitnessQuantity;

    protected:
      //! \brief Evaluator to which all functional calls should be tunneled.
      //! \sa GA::Evaluator, SingleSite::Evaluator, TwoSite::Evaluator,
      //! BandGap::Evaluator, CE::Evaluator, Molecularity::Evaluator
      t_Evaluator &evaluator;
      //! \brief Links quantities and Fitness. \sa Fitness, Objective
      t_Objective &objective; 
      t_Store &store; //!< Stores results!! \sa Store

    public:
      //! records the number of effective calls to the functional
      types::t_unsigned nb_eval; 
      //! records the number of effective calls to the gradient of the functional
      types::t_unsigned nb_grad; 

    public:
      //! \brief Constructor and Initializer
      Base   ( t_Evaluator &_eval, t_Objective &_obj, t_Store &_store )
           : evaluator(_eval), objective(_obj), store(_store), nb_eval(0), nb_grad(0) {};
      //! \brief Copy Constructor 
      Base   ( const Base<t_Evaluator> &_x )
           : evaluator(_x.evaluator), objective(_x.objective),
             store(_x.store), nb_eval(_x.nb_eval), nb_grad(_x.nb_grad) {};
      //! \brief Destructor
      virtual ~Base() {};

    public:
      //! \brief  calls on the functional to evaluate \a _indiv
      //¡ \details All individuals for which the quantities are unknowned
      //! should be evaluated <em>via</em>  this function. The evaluation call
      //! itself is redirected to the Evaluator if and only if the individual is
      //! invalid (meaning its quantities are unknown). The results (e. g. the
      //! quantities) are themselves directed as input to Base::objective. The
      //! result of that is stored into \a _indiv's Fitness object. Finally,
      //! the function also m akes sure that worthy individuals are kept as
      //! results in the object Base::store.
      //! \sa GA::Evaluator::evaluate, SingleSite::Evaluator::evaluate,
      //! TwoSites::Evaluator::evaluate
      virtual t_FitnessQuantity evaluate( t_Individual &_indiv );
      //! \brief calls on the functional to evaluate the (VA) gradients of _indiv
      //! \details Simply links together Objective and the gradients of
      //! Individual::Base::quantities(). Also keeps a tab on the number of
      //! evaluated gradients.
      virtual void evaluate_gradient( t_Individual &_indiv,
                                      t_QuantityGradients& _grad,
                                      t_VA_Type *_i_grad )
      {
        // resets _i_grad
        types::t_unsigned N = _indiv.Object().Container().size();
        t_VA_Type *first = _i_grad;
        std::fill_n( first, N, t_VA_Type(0) );
        // size and value of _grad should be set by evaluator
        evaluator.evaluate_gradient( _grad );
        objective.evaluate_gradient( _indiv.const_quantities(), _grad, _i_grad );
        nb_grad += _grad.size(); 
      }
      //! \brief combines Evaluation::Base::evaluate(t_Individual&) with 
      //! Evaluation::Base::evaluate_gradient().
      virtual t_FitnessQuantity evaluate_with_gradient( t_Individual &_indiv,
                                                        t_QuantityGradients& _grad,
                                                        t_VA_Type *_i_grad );
      //! \brief calls on the functional to evaluate the (VA) gradients of _indiv
      //! \details Simply links together Objective and the gradients of
      //! Individual::Base::quantities(). Also keeps a tab on the number of
      //! evaluated gradients. Only one gradient is evaluated: \a _pos 
      //! \param _indiv Inidividual from which to get gradient
      //! \param _grad Storage for all gradients, as filled by GA::Evaluator derived classes.
      //! \param _pos index of the gradient to evaluate (eg. position in the
      //! presumably vectorial container Traits::Indiv::t_Object::Container)
      virtual t_VA_Type evaluate_one_gradient( t_Individual &_indiv,
                                               t_QuantityGradients& _grad,
                                               types::t_unsigned _pos )
      {
        evaluate( _indiv );
        evaluator.evaluate_one_gradient( _grad, _pos );
        return objective.evaluate_one_gradient( _indiv.const_quantities(), _grad, _pos );
      }

    protected:
      //! \brief Wrapper function which calls Base::evaluate(t_Individual &)
      //! for each individual in \a _pop
      //! \details Additionally, this function first calls Base::init to make
      //! sure that Objective::Base::current_individual and
      //! GA::Evaluator::init are called prior to evaluation. In the _MPI
      //! flavor. it also synchronizes Base::store across all processors.
      virtual void evaluate( t_Population &_pop );
    public:
      //! \brief Always evaluates \a _offspring <em>via</em>
      //! Evaluation::Base::evaluate(t_Population&), conditional evaluation of
      //! \a _pop. 
      //! \details \a _pop is assessed only if the objectives are invalid for
      //! some reason. Probably useless since the introduction of Ranking.
      virtual void operator()( t_Population &_pop, t_Population &_offspring ) 
      { 
        evaluate(_offspring); 
        if ( not objective.is_valid() )
        {
          // if invalid, recomputes whole population
          evaluate( _pop );
          evaluate( _offspring );
        }

        typename t_Population :: const_iterator i_indiv = _offspring.begin();
        typename t_Population :: const_iterator i_end   = _offspring.end();
        Print::out << "New Individual"
                   << ( _offspring.size() > 1 ? "s\n": "   ");
        for(; i_indiv != i_end; ++i_indiv )
        {
          Print::out << *i_indiv << "    Fitness: "
                     << Print::fixed << Print::setw(12) << Print::setprecision(5) << "  "
                     << i_indiv->fitness() << "\n";
        }
        
      }


      //! \brief Correlates Evaluation::Base::objective and
      //! Evaluation::Base::evaluator with \a _indiv
      //! \details This function makes sure that both the Objectives and the
      //! Evaluators know which Individual \a _indiv we are talking about
      void init( t_Individual &_indiv)
      {
        objective.init( _indiv );
        evaluator.init( _indiv );
      }
  };

template< class T_GATRAITS >
  typename Base<T_GATRAITS> :: t_FitnessQuantity
  Base<T_GATRAITS> :: evaluate( t_Individual &_indiv )
  {
    // only computes "expensive" evaluator functionals once!
    if ( _indiv.invalid() ) 
    {
      // fitness AND quantities of _indiv must be valid from here-on
      ++nb_eval;
      evaluator.evaluate();
    }

    _indiv.set_fitness( objective( _indiv.const_quantities() ) );
    store( _indiv );
    return _indiv.fitness();
  }

template< class T_GATRAITS >
  typename Base<T_GATRAITS> :: t_FitnessQuantity
  Base<T_GATRAITS> :: evaluate_with_gradient( t_Individual &_indiv,
                                              t_QuantityGradients& _grad,
                                              t_VA_Type *_i_grad )
  {
    // only computes "expensive" evaluator functionals once!
    nb_grad += _grad.size(); 
    if ( _indiv.invalid() ) 
    {
      // fitness AND quantities of _indiv must be valid from here-on
      ++nb_eval; 
      evaluator.evaluate_with_gradient( _grad );
    }
    else evaluator.evaluate_gradient( _grad );

    _indiv.set_fitness( objective.evaluate_with_gradient( _indiv.const_quantities(),
                                                          _grad, _i_grad ) );
    store( _indiv );
    return _indiv.fitness();
  }

template< class T_GATRAITS >
  void Base<T_GATRAITS> :: evaluate( t_Population &_pop )
  {
    typename t_Population :: iterator i_indiv = _pop.begin();
    typename t_Population :: iterator i_end = _pop.end();
    for(; i_indiv != i_end; ++i_indiv )
    {
      init( *i_indiv );
      evaluate( *i_indiv );
    }

#ifdef _MPI
    store.synchronize();
#endif
  }

  //! \brief Adds History tracking to Evaluation::Base
  //! \details History is managed <em>via</em> the pointer
  //! Evaluation::WithHistory::history to a GA::History object.
  //! History allows us to evaluate any object only once. 
  //! If Evaluation::WithHistory::history is \c NULL then it is ignored.
template< class T_GATRAITS >
  class WithHistory : public Base<T_GATRAITS>
  {
    public:
      //! Should contain all GA definitions \sa Traits::GA
      typedef T_GATRAITS t_GATraits;

    private:
      typedef Base<t_GATraits>                                     t_Base;
      typedef typename t_GATraits::t_Individual                    t_Individual;
      typedef typename t_GATraits :: t_Evaluator                   t_Evaluator;
      typedef typename t_GATraits :: t_QuantityTraits              t_QuantityTraits;
      typedef typename t_QuantityTraits :: t_Quantity              t_Quantity;
      typedef typename t_QuantityTraits :: t_ScalarQuantity        t_ScalarQuantity;
      typedef typename t_GATraits :: t_Population                  t_Population;
      typedef typename Objective :: Types < t_GATraits > :: Vector t_Objective;
      typedef typename t_GATraits :: t_VA_Traits                   t_VATraits;
      typedef typename t_VATraits :: t_QuantityGradients           t_QuantityGradients;
      typedef typename t_VATraits :: t_Type                        t_VA_Type;
      typedef Store :: Base<t_GATraits>                            t_Store;
      typedef GA::History<t_Individual>                            t_History;
      typedef typename t_GATraits :: t_Fitness :: t_ScalarFitness  t_Fitness;
      typedef typename t_Fitness :: t_Quantity                     t_FitnessQuantity;

    protected:
      using t_Base :: evaluator;
      using t_Base :: objective;
      using t_Base :: store;
      using t_Base :: nb_eval;
      using t_Base :: nb_grad;
      //! \brief Keeps tab on previously assessed inviduals
      //! \details If \c NULL, then no history tracking is done, in which case
      //! Base and WithHistory are equivalent.
      t_History *history;


    public:
      using t_Base::evaluate;

    public:
      //! \brief Constructor and Initializer
      WithHistory   ( t_Evaluator &_eval, t_Objective &_obj, t_Store &_store, t_History *_hist )
           : Base<T_GATRAITS>(_eval, _obj, _store), history(_hist) {};
      //! \brief Copy Constructor 
      WithHistory   ( const WithHistory<t_Individual> &_c )
                  : Base<T_GATRAITS>( _c ), history(_c.history ) {}
      //! \brief Destructor
      virtual ~WithHistory() {};

    protected:
      //! \brief  calls on the functional to evaluate \a _indiv
      //! \details Performs exactly as Evaluation::Base::evaluate, except that
      //! if Evaluation::WithHistory::history is non-zero, it does keep tracks of
      //! all individuals evaluated in the past. If \a _indiv is the clone of a
      //! previously evaluated individual, its evaluation in
      //! Evaluation::WithHistory::history is used. 
      //! \sa Evaluation::Base::evaluate 
      virtual t_FitnessQuantity evaluate( t_Individual &_indiv );
      //! \brief combines Evaluation::WithHistory::evaluate(t_Individual&) with 
      //! Evaluation::Base::evaluate_gradient().
      virtual t_FitnessQuantity evaluate_with_gradient( t_Individual &_indiv,
                                                        t_QuantityGradients& _grad,
                                                        t_VA_Type *_i_grad );
    protected:
#ifdef _MPI
      //! \brief adds history synchronization over all processors to
      //! Base::evaluation(t_Population&)
      //! \details Only defined with _MPI.
      virtual void evaluate( t_Population &_pop )
      {
        t_Base::evaluate( _pop );
        if ( history ) history->synchronize();
      }
#endif
  };

template< class T_GATRAITS >
  typename WithHistory<T_GATRAITS> :: t_FitnessQuantity
  WithHistory<T_GATRAITS> :: evaluate( t_Individual &_indiv )
  {
    bool isnot_clone = (history != NULL); // isnot_clone is false if history does not exist
    bool do_evaluate = _indiv.invalid();

    if ( history and  history->clone( _indiv ) )
    {
      isnot_clone = false;
      do_evaluate = false;
    }

    // only evaluates once! from here-on _indiv's fitness should be valid, 
    // and its "quantities" should have been computed
    if ( do_evaluate )
    {
      evaluator.evaluate();
      ++nb_eval;
    }
    _indiv.set_fitness( objective( _indiv.const_quantities() ) ); 

    // isnot_clone is true only if history exists
    // and prior call to history->clone( _indiv ) returned false
    if( isnot_clone ) history->add( _indiv );
    store( _indiv );

    return _indiv.fitness();
  }

template< class T_GATRAITS >
  typename WithHistory<T_GATRAITS> :: t_FitnessQuantity
  WithHistory<T_GATRAITS> :: evaluate_with_gradient( t_Individual &_indiv,
                                                                 t_QuantityGradients& _grad,
                                                                 t_VA_Type *_i_grad )
  {
    bool isnot_clone = (history != NULL); // isnot_clone is false if history does not exist
    bool do_evaluate = _indiv.invalid();

    if ( history and  history->clone( _indiv ) )
    {
      isnot_clone = false;
      do_evaluate = false;
    }

    // only computes "expensive" evaluator functionals once!
    nb_grad += _grad.size(); 
    if ( do_evaluate )
    {
      // fitness AND quantities of _indiv must be valid from here-on
      ++nb_eval; 
      evaluator.evaluate_with_gradient( _grad );
    }
    else evaluator.evaluate_gradient( _grad );

    _indiv.set_fitness( objective.evaluate_with_gradient( _indiv.const_quantities(), _grad, _i_grad ) );
    store( _indiv );
    return _indiv.fitness();
  }

} // namespace Evaluation
/*@}*/

#endif // _RESULTS_H_

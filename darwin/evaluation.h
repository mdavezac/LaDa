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

namespace Evaluation 
{ 
  // abstract base class for results and storage
template< class T_GATRAITS >
  class Base
  {
    public:
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
      t_Evaluator &evaluator;
      t_Objective &objective;
      t_Store &store;

    public:
      types::t_unsigned nb_eval;
      types::t_unsigned nb_grad;

    public:
      Base   ( t_Evaluator &_eval, t_Objective &_obj, t_Store &_store )
           : evaluator(_eval), objective(_obj), store(_store), nb_eval(0), nb_grad(0) {};
      Base   ( const Base<t_Evaluator> &_x )
           : evaluator(_x.evaluator), objective(_x.objective),
             store(_x.store), nb_eval(_x.nb_eval), nb_grad(_x.nb_grad) {};
      virtual ~Base() {};

    public:
      virtual t_FitnessQuantity evaluate( t_Individual &_indiv );
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
      virtual t_FitnessQuantity evaluate_with_gradient( t_Individual &_indiv,
                                                        t_QuantityGradients& _grad,
                                                        t_VA_Type *_i_grad );
      virtual t_VA_Type evaluate_one_gradient( t_Individual &_indiv,
                                               t_QuantityGradients& _grad,
                                               types::t_unsigned _pos )
      {
        evaluate( _indiv );
        evaluator.evaluate_one_gradient( _grad, _pos );
        return objective.evaluate_one_gradient( _indiv.const_quantities(), _grad, _pos );
      }

    protected:
      virtual void evaluate( t_Population &_pop );
    public:
      // this next function makes sure that 
      //   population-dependent fitnesses ( e.g. hamming distance )
      //   moving-traget fitnesses ( e.g. convex hull )
      // are recomputed adequately
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


      // for VA purposes only
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

  // abstract base class for results and storage
template< class T_GATRAITS >
  class WithHistory : public Base<T_GATRAITS>
  {
    public:
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
      t_History *history;

    public:
      WithHistory   ( t_Evaluator &_eval, t_Objective &_obj, t_Store &_store, t_History *_hist )
           : Base<T_GATRAITS>(_eval, _obj, _store), history(_hist) {};
      WithHistory   ( const WithHistory<t_Individual> &_c )
                  : Base<T_GATRAITS>( _c ), history(_c.history ) {}
      virtual ~WithHistory() {};

    protected:
      virtual t_FitnessQuantity evaluate( t_Individual &_indiv );
      virtual t_FitnessQuantity evaluate_with_gradient( t_Individual &_indiv,
                                                        t_QuantityGradients& _grad,
                                                        t_VA_Type *_i_grad );
    protected:
      // we need redefine this only for mpi sync'ing
#ifdef _MPI
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

#endif // _RESULTS_H_

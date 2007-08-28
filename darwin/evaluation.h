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

#include "objective.h"
#include "store.h"
#include "taboos.h"
#include "gatraits.h"

namespace Evaluation 
{ 
  // abstract base class for results and storage
template<class T_EVALUATOR, class T_GATRAITS = Traits::GA<T_EVALUATOR> >
  class Base
  {
    public:
      typedef T_EVALUATOR  t_Evaluator;
      typedef T_GATRAITS   t_GA_Traits;

    protected:
      typedef typename t_GA_Traits :: t_Individual                               t_Individual;
      typedef typename t_GA_Traits :: t_IndivTraits                              t_IndivTraits;
      typedef typename t_GA_Traits :: t_QuantityTraits                           t_QuantityTraits;
      typedef typename t_QuantityTraits :: t_Quantity                            t_Quantity;
      typedef typename t_QuantityTraits :: t_ScalarQuantity                      t_ScalarQuantity;
      typedef typename t_IndivTraits :: t_Population                             t_Population;
      typedef typename Objective :: Types < t_Evaluator, t_GA_Traits > :: Vector t_Objective;
      typedef Store :: Base<t_Evaluator, t_GA_Traits>                            t_Store;
      typedef typename t_IndivTraits :: t_VA_Traits                              t_VATraits;
      typedef typename t_VATraits :: t_QuantityGradients                         t_QuantityGradients;
      typedef typename t_VATraits :: t_Type                                      t_VA_Type;

    protected:
      t_Evaluator &evaluator;
      t_Objective &objective;
      t_Store *store;

    public:
      types::t_unsigned nb_eval;
      types::t_unsigned nb_grad;

    public:
      Base   ( t_Evaluator &_eval, t_Objective &_obj, t_Store *_store )
           : evaluator(_eval), objective(_obj), store(_store), nb_eval(0), nb_grad(0) {};
      Base   ( const Base<t_Evaluator> &_x )
           : evaluator(_x.evaluator), objective(_x.objective),
             store(_x.store), nb_eval(_x.nb_eval), nb_grad(_x.nb_grad) {};
      virtual ~Base() {};

    public:
      virtual t_ScalarQuantity evaluate( t_Individual &_indiv );
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
        objective.evaluate_gradient( _indiv.quantities(), _grad, _i_grad );
        nb_grad += _grad.size(); 
      }
      virtual t_ScalarQuantity evaluate_with_gradient( t_Individual &_indiv,
                                                       t_QuantityGradients& _grad,
                                                       t_VA_Type *_i_grad );
      virtual t_VA_Type evaluate_one_gradient( t_Individual &_indiv,
                                               t_QuantityGradients& _grad,
                                               types::t_unsigned _pos )
      {
        evaluate( _indiv );
        evaluator.evaluate_one_gradient( _grad, _pos );
        return objective.evaluate_one_gradient( _indiv.quantities(), _grad, _pos );
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
        std::cout << "New Individual " << std::endl;
        for(; i_indiv != i_end; ++i_indiv )
        {
          std::string str; str << (const typename t_IndivTraits::t_Object& ) *i_indiv;
          std::cout  << std::setw(12) << std::setprecision(7) << "  "
                     << str << " "
                     << i_indiv->fitness() << std::endl;
        }
        
      }


      // for VA purposes only
      void init( t_Individual &_indiv)
      {
        objective.init( _indiv );
        evaluator.init( _indiv );
      }
  };

template<class T_EVALUATOR, class T_GATRAITS>
  typename Base<T_EVALUATOR, T_GATRAITS> :: t_ScalarQuantity
  Base<T_EVALUATOR,T_GATRAITS> :: evaluate( t_Individual &_indiv )
  {
    // only computes "expensive" evaluator functionals once!
    if ( _indiv.invalid() ) 
    {
      // fitness AND quantities of _indiv must be valid from here-on
      ++nb_eval;
      evaluator.evaluate();
    }

    types::t_real quantity = objective( _indiv.quantities() );
    _indiv.set_fitness( quantity );
    if( store ) (*store)( _indiv );
    return quantity;
  }

template<class T_EVALUATOR, class T_GATRAITS>
  typename Base<T_EVALUATOR, T_GATRAITS> :: t_ScalarQuantity
  Base<T_EVALUATOR,T_GATRAITS> :: evaluate_with_gradient( t_Individual &_indiv,
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

    types::t_real quantity = objective.evaluate_with_gradient( _indiv.quantities(), _grad, _i_grad );
    _indiv.set_fitness( quantity );
    if( store ) (*store)( _indiv );
    return quantity;
  }

template<class T_EVALUATOR, class T_GATRAITS>
  void Base<T_EVALUATOR,T_GATRAITS> :: evaluate( t_Population &_pop )
  {
    typename t_Population :: iterator i_indiv = _pop.begin();
    typename t_Population :: iterator i_end = _pop.end();
    for(; i_indiv != i_end; ++i_indiv )
    {
      init( *i_indiv );
      evaluate( *i_indiv );
    }

#ifdef _MPI
    if( store ) store->synchronize();
#endif
  }

  // abstract base class for results and storage
template<class T_EVALUATOR, class T_GATRAITS = Traits::GA<T_EVALUATOR> >
  class WithHistory : public Base<T_EVALUATOR,T_GATRAITS>
  {
    public:
      typedef T_EVALUATOR  t_Evaluator;
      typedef T_GATRAITS   t_GA_Traits;

    private:
      typedef Base<t_Evaluator, t_GA_Traits> t_Base;
      typedef typename t_GA_Traits :: t_Individual                               t_Individual;
      typedef typename t_GA_Traits :: t_IndivTraits                              t_IndivTraits;
      typedef typename t_GA_Traits :: t_QuantityTraits                           t_QuantityTraits;
      typedef typename t_QuantityTraits :: t_Quantity                            t_Quantity;
      typedef typename t_QuantityTraits :: t_ScalarQuantity                      t_ScalarQuantity;
      typedef typename t_IndivTraits :: t_Population                             t_Population;
      typedef typename Objective :: Types < t_Evaluator, t_GA_Traits > :: Vector t_Objective;
      typedef typename t_IndivTraits :: t_VA_Traits                              t_VATraits;
      typedef typename t_VATraits :: t_QuantityGradients                         t_QuantityGradients;
      typedef typename t_VATraits :: t_Type                                      t_VA_Type;
      typedef Store :: Base<t_Evaluator, t_GA_Traits>                            t_Store;
      typedef darwin::History<t_Individual>                                      t_History;

    protected:
      using t_Base :: evaluator;
      using t_Base :: objective;
      using t_Base :: store;
      using t_Base :: nb_eval;
      using t_Base :: nb_grad;
      t_History *history;

    public:
      WithHistory   ( t_Evaluator &_eval, t_Objective &_obj, t_Store *_store, t_History *_hist )
           : Base<t_Evaluator>(_eval, _obj, _store), history(_hist) {};
      WithHistory   ( const WithHistory<t_Evaluator> &_c )
                  : Base<t_Evaluator>( _c ), history(_c.history ) {}
      virtual ~WithHistory() {};

    protected:
      virtual t_ScalarQuantity evaluate( t_Individual &_indiv );
      virtual t_ScalarQuantity evaluate_with_gradient( t_Individual &_indiv,
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

template<class T_EVALUATOR, class T_GATRAITS>
  typename WithHistory<T_EVALUATOR, T_GATRAITS> :: t_ScalarQuantity
  WithHistory<T_EVALUATOR,T_GATRAITS> :: evaluate( t_Individual &_indiv )
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
    types::t_real quantity = objective( _indiv.quantities() );
    _indiv.set_fitness( quantity );

    // isnot_clone is true only if history exists and prior call to history->clone( _indiv ) returned false
    if( isnot_clone ) history->add( _indiv );
    if ( store ) (*store)( _indiv );

    return quantity;
  }

template<class T_EVALUATOR, class T_GATRAITS>
  typename WithHistory<T_EVALUATOR, T_GATRAITS> :: t_ScalarQuantity
  WithHistory<T_EVALUATOR,T_GATRAITS> :: evaluate_with_gradient( t_Individual &_indiv,
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

    types::t_real quantity = objective.evaluate_with_gradient( _indiv.quantities(), _grad, _i_grad );
    _indiv.set_fitness( quantity );
    if ( store ) (*store)( _indiv );
    return quantity;
  }

} // namespace Evaluation

#endif // _RESULTS_H_

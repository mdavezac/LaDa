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
      typedef typename t_GA_Traits :: t_Individual       t_Individual;
      typedef typename t_GA_Traits :: t_IndivTraits      t_IndivTraits;
      typedef typename t_GA_Traits :: t_QuantityTraits   t_QuantityTraits;
      typedef typename t_QuantityTraits :: t_Quantity    t_Quantity;
      typedef typename t_IndivTraits :: t_Population     t_Population;
      typedef typename Objective :: Base< t_Quantity >   t_Objective;
      typedef Store :: Base<t_Evaluator, t_GA_Traits>    t_Store;

    protected:
      t_Evaluator &evaluator;
      t_Objective &objective;
      t_Store *store;

    public:
      types::t_unsigned nb_eval;

    public:
      Base   ( t_Evaluator &_eval, t_Objective &_obj, t_Store &_store )
           : evaluator(_eval), objective(_obj), store(&_store), nb_eval(0) {};
      Base   ( const Base<t_Evaluator> &_x )
           : evaluator(_x.evaluator), objective(_x.objective),
             store(_x.store), nb_eval(_x.nb_eval) {};
      virtual ~Base() {};

    protected:
      virtual void evaluate( t_Individual &_indiv );
    public:
      virtual void evaluate( t_Population &_pop );
      virtual void operator()( t_Population &_pop, t_Population & ) { evaluate(_pop); }
  };


template<class T_EVALUATOR, class T_GATRAITS>
  void Base<T_EVALUATOR,T_GATRAITS> :: evaluate( t_Individual &_indiv )
  {
    ++nb_eval;
    evaluator.evaluate();
    types::t_real quantity = objective( _indiv.quantities() );
    _indiv.set_fitness( quantity );
    (*store)( _indiv );
  }

template<class T_EVALUATOR, class T_GATRAITS>
  void Base<T_EVALUATOR,T_GATRAITS> :: evaluate( t_Population &_pop )
  {
    typename t_Population :: iterator i_indiv = _pop.begin();
    typename t_Population :: iterator i_end = _pop.end();
    for(; i_indiv != i_end; ++i_indiv )
    {
      if ( not i_indiv->invalid() ) continue;
    
      evaluator.init( *i_indiv );
      evaluate( *i_indiv );
    }

#ifdef _MPI
    store->synchronize();
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
      typedef typename t_GA_Traits :: t_Individual       t_Individual;
      typedef typename t_GA_Traits :: t_IndivTraits      t_IndivTraits;
      typedef typename t_GA_Traits :: t_QuantityTraits   t_QuantityTraits;
      typedef typename t_QuantityTraits :: t_Quantity    t_Quantity;
      typedef typename t_IndivTraits :: t_Population     t_Population;
      typedef typename Objective :: Base< t_Quantity >   t_Objective;
      typedef Store :: Base<t_Evaluator, t_GA_Traits>    t_Store;
      typedef darwin::History<t_Individual> t_History;

    protected:
      using t_Base :: evaluator;
      using t_Base :: objective;
      using t_Base :: store;
      using t_Base :: nb_eval;
      t_History *history;

    public:
      WithHistory   ( t_Evaluator &_eval, t_Objective &_obj, t_Store &_store, t_History *_hist )
           : Base<t_Evaluator>(_eval, _obj, _store), history(_hist) {};
      WithHistory   ( const WithHistory<t_Evaluator> &_c )
                  : Base<t_Evaluator>( _c ), history(_c.history ) {}
      virtual ~WithHistory() {};

    protected:
      virtual void evaluate( t_Individual &_indiv );
    public:
      virtual void evaluate( t_Population &_pop );
  };

template<class T_EVALUATOR, class T_GATRAITS>
  void WithHistory<T_EVALUATOR,T_GATRAITS> :: evaluate( t_Individual &_indiv )
  {
    ++nb_eval;
    evaluator.evaluate();
    types::t_real quantity = objective( _indiv.quantities() );
    _indiv.set_fitness( quantity );

    if( history ) 
      history->add( _indiv );
    (*store)( _indiv );
  }

template<class T_EVALUATOR, class T_GATRAITS>
  void WithHistory<T_EVALUATOR,T_GATRAITS> :: evaluate( t_Population &_pop )
  {
    if ( not history )
    {
      Base<t_Evaluator>::evaluate( _pop );
      return;
    }

    typename t_Population :: iterator i_indiv = _pop.begin();
    typename t_Population :: iterator i_end = _pop.end();
    for(; i_indiv != i_end; ++i_indiv )
    {
      if (    ( not i_indiv->invalid() ) 
           or history->clone( *i_indiv ) ) continue;
    
      evaluator.init( *i_indiv );
      evaluate( *i_indiv );
    }
#ifdef _MPI
    history->synchronize();
#endif
  }

} // namespace Evaluation

#endif // _RESULTS_H_

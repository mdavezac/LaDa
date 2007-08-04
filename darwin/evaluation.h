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

namespace Evaluation 
{ 
  // abstract base class for results and storage
  template<class T_EVALUATOR, class T_POPULATION = eoPop<typename T_EVALUATOR::T_INDIVIDUAL> >
  class Base
  {
    public:
      typedef T_POPULATION t_Population;
      typedef T_EVALUATOR  t_Evaluator;

    protected:
      typedef typename t_Evaluator::t_Individual t_Individual;
      typedef typename Objective :: Base<typename t_Individual::t_Quantity > t_Objective;
      typedef Store :: Base<t_Evaluator> t_Store;

    protected:
      t_Evaluator &evaluator;
      t_Objective &objective;
      t_Store *store;

    public:
      types::t_unsigned nb_eval;

    public:
      Base   ( t_Evaluator &_eval, t_Objective &_obj, t_Store &_store )
           : evaluator(_eval), objective(_obj), store(&_store), nb_eval(0) {};
      virtual ~Base() {};

      void set_evaluator ( t_Evaluator *_e ) { evaluator = _e; }
      
      virtual void evaluate( t_Individual &_indiv );
      virtual void evaluate( t_Population &_pop );
      virtual bool Save( TiXmlElement &_node ) {};
      virtual bool Load( const TiXmlElement &_node ) {};
  };


template<class T_EVALUATOR, class T_POPULATION = eoPop<typename T_EVALUATOR::T_INDIVIDUAL> >
  virtual void Base<T_EVALUATOR, T_POPULATION> :: evaluate( t_Individual &_indiv )
  {
    if ( not _indiv.invalid() )
      return;
    
    ++nb_eval;
    evaluator->evaluate( _indiv );
    types::t_real quantity = objective( _indiv.object.quantities() );
    _indiv.set_fitness( quantity );
    (*store)( _indiv );
  }

template<class T_EVALUATOR, class T_POPULATION = eoPop<typename T_EVALUATOR::T_INDIVIDUAL> >
  virtual void Base<T_EVALUATOR, T_POPULATION> :: evaluate( t_Population &_pop )
  {
    typename t_Population :: iterator i_indiv = _pop.begin();
    typename t_Population :: iterator i_end = _pop.end();
    for(; i_indiv != i_end; ++i_indiv )
    {
      evaluator->init( _indiv );
      evaluate( *i_indiv );
    }

#ifdef _MPI
    store->synchronize();
#endif
  }

  // abstract base class for results and storage
template<class T_EVALUATOR, class T_POPULATION = eoPop<typename T_EVALUATOR::T_INDIVIDUAL> >
  class WithHistory : public Base<T_EVALUATOR, T_POPULATION>
  {
    public:
      typedef T_EVALUATOR  t_Evaluator;
      typedef T_POPULATION t_Population;

    private:
      typedef typename Base<t_Evaluator, t_Population> t_Base;
      typedef typename t_Evaluator :: t_Individual;
      using t_Base :: t_Objective;
      typedef History<t_Individual, std::list<t_Individual> > t_History;
      typedef store :: Base t_Store;

    protected:
      using Base<t_Evaluator, t_Population> :: evaluator;
      using Base<t_Evaluator, t_Population> :: objective;
      using Base<t_Evaluator, t_Population> :: store;
      using Base<t_Evaluator, t_Population> :: nb_eval;
      t_History *history;

    public:
      types::t_unsigned nb_eval;

    public:
      WithHistory   ( t_Evaluator &_eval, t_Objective &_obj, t_Store &_store, t_History *_hist )
           : Base(_eval, _obj, _store), history(_hist) {};
      virtual ~Base() {};

      void set_evaluator ( t_Evaluator *_e ) { evaluator = _e; }
      
      virtual void evaluate( t_Individual &_indiv );
#ifdef _MPI
      virtual void evaluate( t_Population &_pop );
#endif
  };

template<class T_EVALUATOR, class T_POPULATION = eoPop<typename T_EVALUATOR::T_INDIVIDUAL> >
  virtual void WithHistory<T_EVALUATOR, T_POPULATION> :: evaluate( t_Individual &_indiv )
  {
    if ( not _indiv.invalid() )
      return;
    if ( history and history->clone( _indiv ) )
      return;
    
    ++nb_eval;
    evaluator->init( _indiv );
    evaluator->evaluate( _indiv );
    types::t_real quantity = objective( _indiv.object.quantities() );
    _indiv.set_fitness( quantity );

    if( history ) 
      history->add( _indiv );
    (*store)( _indiv );
  }

#ifdef _MPI
template<class T_EVALUATOR, class T_POPULATION = eoPop<typename T_EVALUATOR::T_INDIVIDUAL> >
  virtual void WithHistory<T_EVALUATOR, T_POPULATION> :: evaluate( t_Population &_pop )
  {
    Base<t_Evaluator, t_Population> :: evaluate( _pop );
    if ( history )
      history->synchronize();
  }
#endif

} // namespace Evaluation

#endif // _RESULTS_H_

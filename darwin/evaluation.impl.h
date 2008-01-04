//
//  Version: $Id$
//
#ifndef _EVALUATION_IMPL_H_
#define _EVALUATION_IMPL_H_

namespace Evaluation
{
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

  template< class T_GATRAITS >
    void Base<T_GATRAITS> :: evaluate_gradient( t_Individual &_indiv,
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

  template< class T_GATRAITS >
    typename Base<T_GATRAITS> :: t_VA_Type
      Base<T_GATRAITS> :: evaluate_one_gradient( t_Individual &_indiv,
                                                 t_QuantityGradients& _grad,
                                                 types::t_unsigned _pos )
      {
        ++nb_grad;
//       if( _indiv.invalid() or (not objective.is_valid()) ) evaluate( _indiv );
        evaluate( _indiv );
        evaluator.evaluate_one_gradient( _grad, _pos );
        return objective.evaluate_one_gradient( _indiv.const_quantities(), _grad, _pos );
      }

  template< class T_GATRAITS >
    void Base<T_GATRAITS> :: operator()( t_Population &_pop, t_Population &_offspring ) 
    { 
      evaluate(_offspring); 

      if ( objective.is_valid() ) return;
      
      // if invalid, recomputes whole population
      evaluate( _pop );
      evaluate( _offspring );
    }


  template< class T_GATRAITS >
    inline void Base<T_GATRAITS> :: init( t_Individual &_indiv)
      {
        objective.init( _indiv );
        evaluator.init( _indiv );
      }




#ifdef _MPI
  template< class T_GATRAITS >
    inline void WithHistory<T_GATRAITS> :: evaluate( t_Population &_pop )
    {
      t_Base::evaluate( _pop );
      if ( history ) history->synchronize();
    }
#endif
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

      _indiv.set_fitness( objective.evaluate_with_gradient( _indiv.const_quantities(),
                                                            _grad, _i_grad ) );
      store( _indiv );
      return _indiv.fitness();
    }
}
#endif

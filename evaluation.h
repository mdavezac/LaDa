#ifndef _EVALUATION_H_
#define _EVALUATION_H_

#include <eo/eoEvalFunc.h>
#include <exception>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <string>
using namespace eo;

namespace LaDa 
{
  // Object should be some LaDa::Individual<..,..>
  // Boss should be some LaDa::MotU<..,..>
  template<class Object, class Boss> 
  class Evaluation : public eoEvalFunc<Object> 
  {
    private:
      static Boss *overlord;
      
    public: 
      Evaluation(Boss *_overlord = NULL)
        { overlord = _overlord;}

      void set_overlord( Boss *_overlord )
        { overlord = _overlord; }

      void operator()(Object &_object)
      {
        if ( not overlord )
          throw std::invalid_argument( "pointer overlord unitialised in void Evaluation(Object &_object)" );

        // if quantity does not exist, neither should baseline
        if ( _object.invalid() )
        {
          _object.set_quantity( overlord->evaluate( _object ) );
          _object.set_baseline( overlord->evaluate( _object.get_concentration() ) );
          _object.set_fitness();
          return;
        }

        // everything OK -- validates
        if ( _object.is_baseline_valid() )
          return;

        // baseline should be recomputed
        _object.set_baseline( overlord->evaluate( _object.get_concentration() ) );
        _object.set_fitness();
      }
  };

  
  template<class Object, class Boss> 
    Boss *Evaluation<Object, Boss> :: overlord = NULL;

  template<class Object, class Boss> 
  class Minimization : public eoEvalFunc<Object> 
  {
    private:
      static Boss *overlord;
      
    public: 
      Minimization(Boss *_overlord = NULL)
        { overlord = _overlord;}

      void set_overlord( Boss *_overlord )
        { overlord = _overlord; }

      void operator()(Object &_object)
      {
        if ( not overlord )
          throw std::invalid_argument( "pointer overlord unitialised in void Evaluation(Object &_object)" );

        // if quantity does not exist, neither should baseline
        if ( _object.invalid() )
        {
          _object.set_quantity( overlord->minimize( _object ) );
          _object.set_baseline( overlord->evaluate( _object.get_concentration() ) );
          _object.set_fitness();
          return;
        }

        // everything OK -- validates
        if ( _object.is_baseline_valid() )
          return;

        // baseline should be recomputed
        _object.set_baseline( overlord->evaluate( _object.get_concentration() ) );
        _object.set_fitness();
      }
  };

  template<class Object, class Boss> 
    Boss *Minimization<Object, Boss> :: overlord = NULL;

  template<class Object>
  class EvaluatePop : public eoPopEvalFunc<Object>
  {
    public:
      EvaluatePop(eoEvalFunc<Object> * _eval,
                  eoEvalFunc<Object> *_minimize) : eval(_eval), minimize(_minimize),
                                                   minimize_offsprings(true) {};
    
      void operator()(eoPop<Object> & _parents, eoPop<Object> & _offsprings)
      {
        typename std::vector<Object> :: iterator i_pop = _offsprings.begin();
        typename std::vector<Object> :: iterator i_last = _offsprings.end();

        if ( minimize_offsprings )
          for ( ; i_pop != i_last; ++i_pop )
            minimize->operator()(*i_pop);
        else 
          for ( ; i_pop != i_last; ++i_pop )
            eval->operator()(*i_pop);

        // convex hull has changed => reevaluate
        if ( not  _offsprings.begin()->is_baseline_valid() )
        {
          std::cout << std::endl << "Base line changed" << std::endl; 
          i_pop = _offsprings.begin();
          for ( ; i_pop != i_last; ++i_pop )
            eval->operator()(*i_pop);
          i_pop = _parents.begin();
          i_last = _parents.end();
          for ( ; i_pop != i_last; ++i_pop )
            eval->operator()(*i_pop);
        } 
        i_pop = _offsprings.begin();
        i_last = _offsprings.end();
        std::cout << "New Individuals:" << std::endl; 
        for (int i = 0 ; i_pop != i_last; ++i, ++i_pop )
          std::cout << " Offspring " << i 
                    << " Fitness: " << i_pop->fitness() 
                    << " Quantity: " << i_pop->get_quantity() 
                    << " Baseline: " << i_pop->get_baseline() << std::endl; 
        std::cout << std::endl; 
      }

      void set_minimize_offsprings( bool _b )
        { minimize_offsprings = _b; }
    
    private:
      eoEvalFunc<Object> * eval;
      eoEvalFunc<Object> * minimize;
      bool minimize_offsprings;
  };

  template<class Object>
  class MinimizeBest : public eoPopAlgo<Object>
  {
    public:
      MinimizeBest   (eoEvalFunc<Object> &_minimizer, double _rate = 0.0, unsigned _every = 5)
                   : do_minimize_best(false),
                     how_many( _rate ),
                     minimizer(_minimizer),
                     every(5),
                     nb_calls(0)
      {
       if( _rate > 0 and _rate <= 1 )
         do_minimize_best = true;
       if( _every == 0 ) 
         do_minimize_best = false;
      };

      MinimizeBest   ( const MinimizeBest<Object> &_algo )
                   : do_minimize_best( _algo.do_minimize_best ),
                     how_many( _algo.how_many ),
                     minimizer( _algo.minimizer ),
                     every( _algo.every ),
                     nb_calls( _algo.nb_calls )
                   {}
    
      void operator()(eoPop<Object> & _pop)
      {
        ++nb_calls;
        if ( not do_minimize_best or nb_calls % every )
          return;

        unsigned pSize = _pop.size();
        unsigned nb = how_many( pSize );
        if (nb > pSize )
          return;

        // reorders elements such that the nth best are the nth first
        _pop.nth_element(nb);
        
        typename std::vector<Object> :: iterator i_pop = _pop.begin();
        typename std::vector<Object> :: iterator i_last = _pop.end();
        for (unsigned i = 0; i < nb and i_pop != i_last; ++i, ++i_pop )
        {
          i_pop->invalidate(); // unsets fitness so evaluation is performed
          minimizer( *i_pop );
        }

        if ( not _pop.begin()->is_baseline_valid() )
        {
          i_pop = _pop.begin();
          for ( ; i_pop != i_last; ++i_pop )
            minimizer(*i_pop); // fitness is set, hence won't minimize
        } 

      }

    private:
      bool do_minimize_best;
      eoHowMany how_many;
      eoEvalFunc<Object> & minimizer;
      unsigned every;
      unsigned nb_calls;
  };
  
} // endif LaDa

#endif

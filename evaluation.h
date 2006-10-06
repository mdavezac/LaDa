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


  template<class Object>
  class EvaluatePop : public eoPopEvalFunc<Object>
  {
    public:
      EvaluatePop(eoEvalFunc<Object> & _eval) : eval(_eval) {}
    
      void operator()(eoPop<Object> & _parents, eoPop<Object> & _offsprings)
      {
        typename std::vector<Object> :: iterator i_pop = _offsprings.begin();
        typename std::vector<Object> :: iterator i_last = _offsprings.end();
        for ( ; i_pop != i_last; ++i_pop )
          eval(*i_pop);
        // convex hull has changed => reevaluate
        if ( not  _offsprings.begin()->is_baseline_valid() )
        {
          std::cout << std::endl << "Base line changed" << std::endl; 
          i_pop = _offsprings.begin();
          for ( ; i_pop != i_last; ++i_pop )
            eval(*i_pop);
          i_pop = _parents.begin();
          i_last = _parents.end();
          for ( ; i_pop != i_last; ++i_pop )
            eval(*i_pop);
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
    
    private:
      eoEvalFunc<Object> & eval;
  };
  
} // endif LaDa

#endif

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
          i_pop = _offsprings.begin();
          for ( ; i_pop != i_last; ++i_pop )
            eval(*i_pop);
          i_pop = _parents.begin();
          i_last = _parents.end();
          for ( ; i_pop != i_last; ++i_pop )
            eval(*i_pop);
        } 
       //  int i=0, j=0;
       //try
       //{ 
       //  for(; i < _offsprings.size(); ++i)
       //    _offsprings[i].fitness();
       //  for(; j < _parents.size(); ++j)
       //    _parents[j].fitness();
       //}
       //catch (std::exception &e )
       //{
       //  std::cerr << "here" << std::endl;
       //  std::ostringstream s( e.what() );
       //  s << " at i=" << i << " and j=" << j <<" ";
       //  throw std::runtime_error(s.str());
       //}
      }
    
    private:
      eoEvalFunc<Object> & eval;
  };
  
} // endif LaDa

#endif

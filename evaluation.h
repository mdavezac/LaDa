#ifndef _EVALUATION_H_
#define _EVALUATION_H_

#include <eo/eoEvalFunc.h>
#include <exception>
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
        if ( not _object.invalid() )
          return;

        if ( not overlord )
          throw std::invalid_argument( "pointer overlord unitialised in void Evaluation(Object &_object)" );

        if ( not _object.is_quantity_valid() )
          _object.set_quantity( overlord->evaluate( _object ) );

        if ( not _object.is_baseline_valid() )
        {
          _object.set_baseline( overlord->evaluate( _object.get_concentration() ) );
          _object.set_fitness();
          return;
        }

        if ( _object.invalid() )
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
        apply<Object>(eval, _offsprings);
        // convex hull has changed => reevaluate
        if ( not  _offsprings.begin()->is_baseline_valid() )
        {
          apply<Object>(eval, _offsprings);
          apply<Object>(eval, _parents);
        }
      }
    
    private:
      eoEvalFunc<Object> & eval;
  };
  
} // endif LaDa

#endif

#ifndef _EVALUATION_H_
#define _EVALUATION_H_

#include <eo/eoEvalFunc.h>
#include <exception>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <string>
using namespace eo;

#include "operators.h" 
#include <algorithm>
#include <string>

#include <opt/types.h>
using namespace types;
#include "eotypes.h"

namespace LaDa 
{
  // t_Object should be some LaDa::Individual<..,..>
  // Boss should be some LaDa::MotU<..,..>
  template<class t_Object, class Boss> 
  class Evaluation : public eoEvalFunc<t_Object> 
  {
    private:
      Boss &overlord;
      
    public: 
      Evaluation(Boss &_overlord) : overlord( _overlord ) {}

      void set_overlord( Boss &_overlord )
        { overlord = _overlord; }

      void operator()(t_Object &_object)
      {
        // if quantity does not exist, neither should baseline
        if ( _object.invalid() )
        {
          _object.set_quantity( overlord.evaluate( _object ) );
          _object.set_baseline( overlord.evaluate( _object.get_concentration() ) );
          _object.set_fitness();
          return;
        }

        // everything t_Object -- validates
        if ( _object.is_baseline_valid() )
          return;

        // baseline should be recomputed
        _object.set_baseline( overlord.evaluate( _object.get_concentration() ) );
        _object.set_fitness();
      }
  };

  
  template<class t_Object>
  class EvaluatePop : public eoPopEvalFunc<t_Object>
  {
    public:
      EvaluatePop(eoEvalFunc<t_Object> &_eval ) : eval(_eval) {}
    
      void operator()(eoPop<t_Object> & _parents, eoPop<t_Object> & _offsprings)
      {
        typename std::vector<t_Object> :: iterator i_pop = _offsprings.begin();
        typename std::vector<t_Object> :: iterator i_last = _offsprings.end();

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
        for (t_int i = 0 ; i_pop != i_last; ++i, ++i_pop )
          std::cout << " Offspring " << i 
                    << " Fitness: " << i_pop->fitness() 
                    << " Quantity: " << i_pop->get_quantity() 
                    << " Baseline: " << i_pop->get_baseline() << std::endl; 
        std::cout << std::endl; 
      }

    private:
      eoEvalFunc<t_Object> &eval;
  };

  template<class t_Object, class t_Call_Back>
  class Extra_PopAlgo : public eoPopAlgo<t_Object>
  {
      bool do_minimize_best;
      t_Call_Back &call_back;
      eoHowMany how_many;
      t_unsigned every;
      t_unsigned nb_calls;
      Evaluation<t_Object, t_Call_Back> *eval;
      eoMonOp<t_Object> & op;

    public:
      Extra_PopAlgo  (eoMonOp<t_Object> &_op, t_Call_Back &_call_back,
                      t_real _rate = 0.0, t_unsigned _every = 5)
                   : do_minimize_best(false),
                     call_back(_call_back), 
                     how_many( _rate ),
                     every(_every),
                     nb_calls(0),
                     op(_op)
      {
        eval = new Evaluation<t_Object, t_Call_Back>(call_back);
        if( _rate > 0 and _rate <= 1 )
          do_minimize_best = true;
        if( _every == 0 ) 
          do_minimize_best = false;
      };

      Extra_PopAlgo  ( const Extra_PopAlgo<t_Object, t_Call_Back> &_algo )
                   : do_minimize_best( _algo.do_minimize_best ),
                     call_back(_algo.call_back), 
                     how_many( _algo.how_many ),
                     op( _algo.op ),
                     every( _algo.every ),
                     nb_calls( _algo.nb_calls ),
                     eval(_algo.eval),
                     op(_algo.op)
                   {}
      ~Extra_PopAlgo () { delete eval; }
    
      void operator()(eoPop<t_Object> & _pop)
      {
        ++nb_calls;
        if ( not do_minimize_best or nb_calls % every )
          return;

        t_unsigned pSize = _pop.size();
        eotypes::t_unsigned nb = how_many( (eotypes::t_unsigned)pSize );
        if (nb > pSize )
          return;

        {
          std::string str(" ExtraPopAlgo launched ");
          call_back.print_xmgrace( str );
        }


        // reorders elements such that the nth best are the nth first
        _pop.nth_element(nb);
        
        typename std::vector<t_Object> :: iterator i_pop = _pop.begin();
        typename std::vector<t_Object> :: iterator i_last = _pop.end();
        for (t_unsigned i = 0; i < nb and i_pop != i_last; ++i, ++i_pop )
        {
          i_pop->invalidate();
          op( *i_pop );
          (*eval)(*i_pop);
        }

        if ( not _pop.begin()->is_baseline_valid() )
        {
          i_pop = _pop.begin();
          for (; i_pop != i_last; ++i_pop )
            (*eval)(*i_pop);
        }

      }

  };
  
} // endif LaDa

#endif

//
//  Version: $Id$
//
#ifndef _DARWIN_OPERATORS_H_
#define _DARWIN_OPERATORS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <exception>
#include <iostream>
#include <vector>
#include <list>
#include <utility>

#include <eo/eoPop.h>
#include <eo/eoOp.h>
#include <eo/eoOpContainer.h>
#include <eo/utils/eoRNG.h>



#include "opt/types.h"
#include "gencount.h"

namespace GA 
{
  template<class T_GATRAITS>
  class SequentialOp : public eoOpContainer<typename T_GATRAITS :: t_Individual >
  {
    protected:
      typedef T_GATRAITS t_GATraits;
    protected:
      typedef typename t_GATraits :: t_Individual t_Individual;
      using eoOpContainer< t_Individual >::ops;
      using eoOpContainer< t_Individual >::rates;

    public:
      SequentialOp() {}

      virtual void apply(eoPopulator<t_Individual> &_populator)
      {
        typename std::vector< types::t_real > :: const_iterator i_rate = rates.begin();
        typename std::vector< eoGenOp<t_Individual>* > :: iterator i_op = ops.begin();
        typename std::vector< eoGenOp<t_Individual>* > :: iterator i_end = ops.end();
        for (; i_op != i_end; ++i_op, ++i_rate )
        {
          if ( eo::rng.flip( *i_rate ) )
            (*i_op)->operator()( _populator );
        }
      }
    
      virtual std::string className() const {return "GA::SequentialOp";}

  };

  template<class T_GATRAITS>
  class ProportionalOp : public eoOpContainer<typename T_GATRAITS :: t_Individual >
  {
    protected:
      typedef T_GATRAITS t_GATraits;
    protected:
      typedef typename t_GATraits :: t_Individual t_Individual;
      using eoOpContainer< t_Individual >::ops;
      using eoOpContainer< t_Individual >::rates;

    public:
      virtual void apply(eoPopulator<t_Individual> &_populator)
      {
        types::t_unsigned i = static_cast<types::t_unsigned>(eo::rng.roulette_wheel(rates));
        return (*ops[i])( _populator );
      }

      virtual std::string className() const {return "GA::ProportionalOp";}
  };

  template<class T_GATRAITS>
  class TriggeredOp : public eoGenOp<typename T_GATRAITS :: t_Individual >
  {
    public:
      typedef T_GATRAITS t_GATraits;
    protected:
      typedef typename t_GATraits :: t_Individual t_Individual;

    protected:
      bool is_triggered;
      eoGenOp<t_Individual> *op;
      ;
    public:
      TriggeredOp  ( eoOp<t_Individual> &_op,
                     eoFunctorStore &_store,
                     bool _t = false  )
                  : is_triggered(_t)
        { op = &wrap_op<t_Individual>( _op, _store ); }
      virtual ~TriggeredOp() {};
    
      virtual types::t_unsigned max_production()
        { return op->max_production(); }

      void trigger( bool _trigger = false )
      {
        is_triggered = _trigger;
      }
      virtual std::string className() const {return "GA :: TriggeredOps";}

      virtual void apply( eoPopulator<t_Individual> &_indiv ) 
      {
        if ( is_triggered )
          (*op)(_indiv);
      }

  };

  template<class T_GATRAITS>
  class PeriodicOp : public eoGenOp<typename T_GATRAITS :: t_Individual >
  {
    public:
      typedef T_GATRAITS t_GATraits;
    protected:
      typedef typename t_GATraits :: t_Individual t_Individual;
      
    protected:
      types::t_unsigned period;
      GenCount &age;
      eoGenOp<t_Individual> *op;
    public:
      PeriodicOp   ( eoOp<t_Individual> &_op,
                     types::t_unsigned _period, 
                     GenCount &_age,
                     eoFunctorStore &_store)
                 : period(_period), age(_age)
        { op = &wrap_op<t_Individual>( _op, _store ); }
      virtual ~PeriodicOp() {};
    
      virtual types::t_unsigned max_production()
        { return op->max_production(); }

      virtual std::string className() const {return "GA :: PeriodicOp";}

      virtual void apply( eoPopulator<t_Individual> &_indiv ) 
      {
        if ( age() % period == 0 )
          (*op)(_indiv);
      }
  };


  template<class T_GATRAITS>
  class AgeOp : public eoMonOp<typename T_GATRAITS :: t_Individual >
  {
    public:
      typedef T_GATRAITS t_GATraits;
    protected:
      typedef typename t_GATraits :: t_Individual t_Individual;
      
    protected:
      GenCount &age;

    public:
      AgeOp   ( GenCount &_age )
            : age(_age) {}
      AgeOp   ( const AgeOp<t_Individual> &_t )
            : age(_t.age) {}

      // tries to create an untaboo indiv on applying _op
      // after max tries, creates a random untaboo indiv
      virtual bool operator()( t_Individual &_indiv )
      {
        _indiv.set_age( age() );
        return false;
      }
  };
} // endif GA

#endif

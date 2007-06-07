#ifndef _OPERATORS_H_
#define _OPERATORS_H_
//   defines a mutation and a crossover operator
//   t_Object template should be some LaDa::Individual<..,..>
//   also defines a randomization operator

#include <exception>
#include <iostream>

#include <eo/eoPop.h>
#include <eo/eoOp.h>
#include <eo/eoOpContainer.h>
#include <eo/utils/eoRNG.h>

#include <opt/opt_minimize.h>

#include<vector>
#include<list>
#include<utility>

#include <opt/types.h>
using namespace types;
#include "eotypes.h"
#include "gencount.h"

namespace LaDa 
{
  template<class t_Object> 
  class Crossover : public eoBinOp<t_Object> 
  {
    private:
      t_real probability;

    public:
      Crossover( t_real _c = 0.5 ) : probability(_c) {};

      void set_probability(t_real &_c)
        { probability = _c; }

      virtual std::string className() const { return "LaDa::Crossover"; }

      bool operator() (t_Object &object1, const t_Object &object2) 
      {
        typename t_Object :: iterator i_object1, i_last;
        typename t_Object :: const_iterator i_object2;

        i_object1 = object1.begin();
        i_object2 = object2.begin();
        i_last = object1.end();
        for( ; i_object1 != i_last; ++i_object1, ++i_object2 )
          if ( rng.uniform() < probability ) 
            *i_object1 = *i_object2;

        if ( t_Object :: is_using_phenotype )
          object1.set_phenotype_to_genotype();
        object1.invalidate();
        
        return true;
      }

  }; // class Crossover : public eoGenOp<t_Object>

// template<class t_Call_Back, class t_Object = typename t_Call_Back :: t_Object >
// class Krossover : public eoBinOp<t_Object> 
// {
//   private:
//     t_Call_Back &call_back;
//
//   public:
//     Krossover( t_Call_Back &_call_back) : call_back(_call_back) {};
//
//     virtual std::string className() const { return "LaDa::Krossover"; }
//
//     bool operator() (t_Object &_object1, const t_Object &_object2) 
//       {  return call_back.Krossover( _object1, _object2 ); }
//
// }; // class Krossover : public eoGenOp<t_Object>
  
  template<class t_Object> 
  class Mutation : public eoMonOp<t_Object> 
  {
    private:
      t_real probability;

    public:
      Mutation( t_real _c = 0.0 ) : probability(_c) {};

      virtual std::string className() const { return "LaDa::Mutation"; }

      void set_probability(t_real &_c)
        { probability = _c; }

      bool operator() (t_Object &object) 
      {
        typename t_Object :: iterator i_object, i_last;
        bool mutated = false;

        i_object = object.begin();
        i_last = object.end();
        for( ; i_object != i_last; ++i_object )
          if ( rng.uniform() < probability ) 
          {
            *i_object *= -1.0;
            mutated = true;
          }

        if ( mutated )
          object.invalidate();

        if ( t_Object :: is_using_phenotype )
          object.set_phenotype_to_genotype();
        
        return mutated;
      }

  }; // class Mutation<t_Object> : public eoMonOp<t_Object> 

  template<class t_Object> 
  class UtterRandom : public eoMonOp<t_Object> 
  {
    public:
      UtterRandom(){};

      virtual std::string className() const { return "LaDa::UtterRandom"; }

      bool operator() (t_Object &object)
      {
        typename t_Object :: iterator i_object, i_last;

        i_object = object.begin();
        i_last = object.end();
        for( ; i_object != i_last; ++i_object )
          *i_object = ( rng.flip() ) ? -1.0 : 1.0;

        object.invalidate();
        
        return true;
      }
    // void print_out( std::ostream &_str)
    //   {  std::cout << "UtterRandom " << probability << " "; }
  }; // class Mutation<t_Object> : public eoMonOp<t_Object> 

  template<class t_Call_Back, class t_Object = typename t_Call_Back :: t_Object > 
  class EvaluateOp : public eoMonOp<t_Object> 
  {
    private:
      t_Call_Back &call_back;

    public:
      EvaluateOp ( t_Call_Back &_call_back ) : call_back(_call_back) {};

      virtual std::string className() const { return "LaDa::EvaluateOp"; }

      bool operator() (t_Object &_object) 
      {  
        call_back.evaluate( _object ); 
        return false;
      }

    // void print_out( std::ostream &_str)
    //   {  std::cout << "Mutation, v=" << probability << " "; }
  }; // class EvaluateOp<t_Object, t_Call_Back> : public eoMonOp<t_Object> 

  template<class EO_OBJECT, class CALL_BACK> 
  class MinimizationOp : public eoMonOp<EO_OBJECT> 
  {
    private:
      t_unsigned minimizer_nb;
      CALL_BACK &call_back;

    public:
      MinimizationOp   ( const MinimizationOp<EO_OBJECT, CALL_BACK> &_minop )
                     : minimizer_nb(_minop.minimizer_nb),
                       call_back( _minop.call_back ) {};
      MinimizationOp   ( t_unsigned _nb, CALL_BACK &_call_back )
                     : minimizer_nb(_nb),
                       call_back( _call_back ) {};
      virtual ~MinimizationOp() {}

  //   void print_out( std::ostream &_str)
  //     {  std::cout << "Minimization "; call_back->print_out_minimizer(minimzer_nb);  }


      virtual std::string className() const { return "LaDa::MinimizerOp"; }

      bool operator() (EO_OBJECT &_object) 
      {
        call_back.minimize( _object, minimizer_nb );
       
        _object.invalidate(); 

        return true;
      }
  }; // class MinimizationOp : public eoMonOp<EO_OBJECT> 

  template<class t_Object>
  class SequentialOp : public eoOpContainer<t_Object>
  {
    protected:
      using eoOpContainer< t_Object >::ops;
      using eoOpContainer< t_Object >::rates;

    public:
      SequentialOp() {}

      virtual void apply(eoPopulator<t_Object> &_populator)
      {
        typename std::vector< eotypes::t_real > :: const_iterator i_rate = rates.begin();
        typename std::vector< eoGenOp<t_Object>* > :: iterator i_op = ops.begin();
        typename std::vector< eoGenOp<t_Object>* > :: iterator i_end = ops.end();
        for (; i_op != i_end; ++i_op, ++i_rate )
        {
          if ( eo::rng.flip( *i_rate ) )
            (*i_op)->operator()( _populator );
        }
      }
    
      virtual std::string className() const {return "LaDa::SequentialOp";}

  };

  template<class t_Object>
  class ProportionalOp : public eoOpContainer<t_Object>
  {
    protected:
      using eoOpContainer< t_Object >::ops;
      using eoOpContainer< t_Object >::rates;

    public:
      virtual void apply(eoPopulator<t_Object> &_populator)
      {
        t_unsigned i = static_cast<t_unsigned>(eo::rng.roulette_wheel(rates));
        return (*ops[i])( _populator );
      }

      virtual std::string className() const {return "LaDa::ProportionalOp";}
  };

  template<class t_Object>
  class TriggeredOp : public eoGenOp<t_Object>
  {
    protected:
      bool is_triggered;
      eoGenOp<t_Object> *op;
      ;
    public:
      TriggeredOp  ( eoOp<t_Object> &_op,
                     eoFunctorStore &_store,
                     bool _t = false  )
                  : is_triggered(_t)
        { op = &wrap_op<t_Object>( _op, _store ); }
      virtual ~TriggeredOp() {};
    
      virtual eotypes::t_unsigned max_production()
        { return op->max_production(); }

      void trigger( bool _trigger = false )
      {
        is_triggered = _trigger;
      }
      virtual std::string className() const {return "LaDa :: TriggeredOps";}

      virtual void apply( eoPopulator<t_Object> &_object ) 
      {
        if ( is_triggered )
          (*op)(_object);
      }

 //    void print_out( std::ostream &_str)
 //    { 
 //      std::cout << "Triggered "; 
 //    }
  };

  template<class t_Object>
  class PeriodicOp : public eoGenOp<t_Object>
  {
    protected:
      t_unsigned period;
      GenCount &age;
      eoGenOp<t_Object> *op;
    public:
      PeriodicOp   ( eoOp<t_Object> &_op,
                     t_unsigned _period, 
                     GenCount &_age,
                     eoFunctorStore &_store)
                 : period(_period), age(_age)
        { op = &wrap_op<t_Object>( _op, _store ); }
      virtual ~PeriodicOp() {};
    
      virtual eotypes::t_unsigned max_production()
        { return op->max_production(); }

      virtual std::string className() const {return "LaDa :: TriggeredOps";}

      virtual void apply( eoPopulator<t_Object> &_object ) 
      {
        if ( age() % period == 0 )
          (*op)(_object);
      }
  };


  template<class t_Object>
  class AgeOp : public eoMonOp<t_Object>
  {
    protected:
      GenCount &age;

    public:
      AgeOp   ( GenCount &_age )
            : age(_age) {}
      AgeOp   ( const AgeOp<t_Object> &_t )
            : age(_t.age) {}

      // tries to create an untaboo object on applying _op
      // after max tries, creates a random untaboo object
      virtual bool operator()( t_Object &_object )
      {
        _object.set_age( age() );
        return false;
      }
  };
} // endif LaDa

#endif

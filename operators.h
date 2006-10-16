#ifndef _OPERATORS_H_
#define _OPERATORS_H_
//   defines a mutation and a crossover operator
//   Object template should be some LaDa::Individual<..,..>
//   also defines a randomization operator

#include <exception>
#include <iostream>

#include <eo/eoPop.h>
#include <eo/utils/eoRNG.h>

#include <opt/opt_minimize.h>

#include<vector>
#include<list>
#include<utility>

namespace LaDa 
{
  template<class Object> 
  class Crossover : public eoBinOp<Object> 
  {
    private:
      double probability;

    public:
      Crossover( double _c = 0.5 ) : probability(_c) {};

      void set_probability(double &_c)
        { probability = _c; }

      virtual std::string className() const { return "LaDa::Crossover"; }

      bool operator() (Object &object1, const Object &object2) 
      {
        typename Object::CONTAINER_ITERATOR i_object1, i_last;
        typename Object::CONST_CONTAINER_ITERATOR i_object2;

        i_object1 = object1.begin();
        i_object2 = object2.begin();
        i_last = object1.end();
        for( ; i_object1 != i_last; ++i_object1, ++i_object2 )
          if ( rng.uniform() < probability ) 
            *i_object1 = *i_object2;

        object1.invalidate();
        
        return true;
      }
  }; // class Crossover : public eoGenOp<Object>

  template<class Object> 
  class Mutation : public eoMonOp<Object> 
  {
    private:
      double probability;

    public:
      Mutation( double _c = 0.0 ) : probability(_c) {};

      virtual std::string className() const { return "LaDa::Mutation"; }

      void set_probability(double &_c)
        { probability = _c; }

      bool operator() (Object &object) 
      {
        typename Object::CONTAINER_ITERATOR i_object, i_last;
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
        
        return mutated;
      }
  }; // class Mutation<Object> : public eoMonOp<Object> 

  template<class Object> 
  class UtterRandom : public eoMonOp<Object> 
  {
    public:
      UtterRandom(){};

      virtual std::string className() const { return "LaDa::UtterRandom"; }

      bool operator() (Object &object)
      {
        typename Object::CONTAINER_ITERATOR i_object, i_last;

        i_object = object.begin();
        i_last = object.end();
        for( ; i_object != i_last; ++i_object )
          *i_object = ( rng.flip() ) ? -1.0 : 1.0;

        object.invalidate();
        
        return true;
      }
  }; // class Mutation<Object> : public eoMonOp<Object> 
  
  template<class EO_OBJECT, class CALL_BACK> 
  class MinimizationOp : public eoMonOp<EO_OBJECT> 
  {
    private:
      unsigned minimizer_nb;
      CALL_BACK &call_back;

    public:
      MinimizationOp   ( const MinimizationOp<EO_OBJECT, CALL_BACK> &_minop )
                     : minimizer_nb(_minop.minimizer_nb),
                       call_back( _minop.call_back ) {};
      MinimizationOp   ( unsigned _nb, CALL_BACK &_call_back )
                     : minimizer_nb(_nb),
                       call_back( _call_back ) {};
      virtual ~MinimizationOp() {}


    protected:
      virtual std::string className() const { return "LaDa::MinimizerOp"; }

      bool operator() (EO_OBJECT &_object) 
      {
        call_back.minimize( _object, minimizer_nb );
       
        _object.invalidate(); 

        return true;
      }
  }; // class MinimizationOp : public eoMonOp<EO_OBJECT> 

  template<class t_Object>
  class SequentialMonOp : public eoMonOp<t_Object>
  {
      typedef std::pair< eoMonOp<t_Object>*, double > t_pair;
      typedef std::list< t_pair > t_pair_list;
      t_pair_list operators;

    public:
      SequentialMonOp() {}

      SequentialMonOp   ( const SequentialMonOp<t_Object> &_algo ) 
                      : operators(_algo.operators ) {};
    
      virtual bool operator()(t_Object &_object)
      {
        typename t_pair_list :: iterator i_op = operators.begin();
        typename t_pair_list :: iterator i_end = operators.end();
        for (; i_op != i_end; ++i_op )
        {
          if ( eo::rng.flip( i_op->second ) )
            (i_op->first)->operator()( _object );
        }

        return false;
      }

      void add( eoMonOp<t_Object> *_op, double &_rate)
      { 
        operators.push_back( t_pair(_op, _rate) ); 
      }
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
                     bool _t = false, eoFunctorStore &_store )
                  : is_triggered(_t)
        { op = &wrap_op<t_Object>( _op, _store ); }
      virtual ~TriggeredOp() {};
    
      virtual unsigned max_production()
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
  };

  template<class t_Object>
  class PeriodicOp : public eoGenOp<t_Object>
  {
    protected:
      unsigned period;
      eoIncrementorParam<unsigned> &age;
      eoGenOp<t_Object> *op;
    public:
      PeriodicOp   ( eoOp<t_Object> &_op,
                     unsigned _period, 
                     eoIncrementorParam<unsigned> &_age,
                     eoFunctorStore &_store)
                 : period(_period), age(_age)
      
        { op = &wrap_op<t_Object>( _op, _store ); }
      virtual ~PeriodicOp() {};
    
      virtual unsigned max_production()
        { return op->max_production(); }

      virtual std::string className() const {return "LaDa :: TriggeredOps";}

      virtual void apply( eoPopulator<t_Object> &_object ) 
      {
        if ( age.value() % period == 0 )
          (*op)(_object);
      }
  };

  template<class t_Object>
  class ProportionalMonOp : public eoMonOp<t_Object>
  {
      std::vector< eoMonOp<t_Object>* > operators;
      std::vector< double > rates;

    public:
      ProportionalMonOp() {}

      ProportionalMonOp   ( const ProportionalMonOp<t_Object> &_algo )
                        : operators(_algo.operators ),
                          rates(_algo.rates) {}
    
      virtual bool operator()(t_Object &_object)
      {
        unsigned i = eo::rng.roulette_wheel(rates);
        return (operators[i])->operator()( _object );
      }

      void add( eoMonOp<t_Object> *_op, double &_rate)
      { 
        operators.push_back( _op ); 
        rates.push_back( _rate );
      }
  };

  template<class t_Object>
  class AgeOp : public eoMonOp<t_Object>
  {
    protected:
      eoIncrementorParam<unsigned> &age;

    public:
      AgeOp   ( eoIncrementorParam<unsigned> &_age )
            : age(_age) {}
      AgeOp   ( const AgeOp<t_Object> &_t )
            : age(_t.age) {}

      // tries to create an untaboo object on applying _op
      // after max tries, creates a random untaboo object
      virtual bool operator()( t_Object &_object )
      {
        _object.set_age( age.value() );
        return false;
      }
  };
} // endif LaDa

#endif

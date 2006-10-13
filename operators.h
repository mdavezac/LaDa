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
  class Sequential : public eoMonOp<t_Object>
  {
      std::vector< eoMonOp<t_Object>* > operators;
      std::vector< double > rates;

    public:
      Sequential() {}

      Sequential ( const Sequential<t_Object> &_algo ) : operators(_algo.operators ) {}
    
      virtual bool operator()(t_Object &_pop)
      {
        typename std::vector< double > :: iterator i_rate = rates.begin();
        typename std::vector< eoMonOp<t_Object>* > :: iterator i_op = operators.begin();
        typename std::vector< eoMonOp<t_Object>* > :: iterator i_end = operators.end();
        for (; i_op != i_end; ++i_op, ++i_rate )
        {
          if ( eo::rng.flip( *i_rate ) )
            (*i_op)->operator()( _pop );
        }

        return false;
      }

      void add( eoMonOp<t_Object> *_op, double &_rate)
      { 
        operators.push_back( _op ); 
        rates.push_back( _rate );
      }
  };

  template<class t_Object>
  class Proportional : public eoMonOp<t_Object>
  {
      std::vector< eoMonOp<t_Object>* > operators;
      std::vector< double > rates;

    public:
      Proportional() {}

      Proportional ( const Proportional<t_Object> &_algo ) : operators(_algo.operators ) {}
    
      virtual bool operator()(t_Object &_pop)
      {
        unsigned i = eo::rng.roulette_wheel(rates);
        return (operators[i])->operator()( _pop );
      }

      void add( eoMonOp<t_Object> *_op, double &_rate)
      { 
        operators.push_back( _op ); 
        rates.push_back( _rate );
      }
  };
} // endif LaDa

#endif

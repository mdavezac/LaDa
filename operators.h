#ifndef _OPERATORS_H_
#define _OPERATORS_H_
//   defines a mutation and a crossover operator
//   Object template should be some LaDa::Individual<..,..>
//   also defines a randomization operator

#include <exception>
#include <iostream>

#include <eo/eoPop.h>
#include <eo/utils/eoRNG.h>

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

        object1.invalidate_quantity();
        
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
          object.invalidate_quantity();
        
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

        
        return true;
      }
  }; // class Mutation<Object> : public eoMonOp<Object> 
  
} // endif LaDa

#endif

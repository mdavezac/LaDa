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
  
  template<class EO_OBJECT, class OPT_OBJECT> 
  class MinimizationOp : public eoMonOp<EO_OBJECT> 
  {
    public:
      // following should be exactly quivalent to MotU :: GA :: ...
      const static unsigned WANG_MINIMIZER; 
      const static unsigned PHYSICAL_MINIMIZER;
      const static unsigned LINEAR_MINIMIZER;
      const static unsigned SA_MINIMIZER; 

    private:
      unsigned n, type;
      OPT_OBJECT *object;
      opt::Minimize_Base<OPT_OBJECT> *minimizer;

    public:
      MinimizationOp   ( const MinimizationOp<EO_OBJECT, OPT_OBJECT> &_minop )
                     : n(_minop.n), type(_minop.type ), object(_minop.object )
        { create_minimizer(); }
      MinimizationOp   ( unsigned _n, unsigned _type, OPT_OBJECT *_fitness )
                     : n(_n), type(_type), object( _fitness )
        { create_minimizer(); };
      virtual ~MinimizationOp()
      {
        if ( minimizer )
          delete minimizer;
        minimizer = NULL;
      }


    protected:
      void create_minimizer()
      {
        switch( type )
        {
          case MinimizationOp<EO_OBJECT, OPT_OBJECT> :: WANG_MINIMIZER: 
            minimizer = new opt::Minimize_Wang<OPT_OBJECT>(object);
            break;
          case MinimizationOp<EO_OBJECT, OPT_OBJECT> :: PHYSICAL_MINIMIZER: 
            minimizer = new opt::Minimize_Ssquared<OPT_OBJECT>(object);
            break;
          case MinimizationOp<EO_OBJECT, OPT_OBJECT> :: SA_MINIMIZER: 
            minimizer = new opt::Minimize_Linear<OPT_OBJECT>(object);
            static_cast< opt::Minimize_Linear<OPT_OBJECT>* >(minimizer)->simulated_annealing = true;
            static_cast< opt::Minimize_Linear<OPT_OBJECT>* >(minimizer)->max_calls = n;
            break;
          case MinimizationOp<EO_OBJECT, OPT_OBJECT> :: LINEAR_MINIMIZER: 
            minimizer = new opt::Minimize_Linear<OPT_OBJECT>(object);
            static_cast< opt::Minimize_Linear<OPT_OBJECT>* >(minimizer)->simulated_annealing = false;
            static_cast< opt::Minimize_Linear<OPT_OBJECT>* >(minimizer)->max_calls = n;
            break;
          default:
            std::cerr << "Unknown minimizer type in LaDa :: MinimizerOp" << std::endl;
            exit(-1);
            break;
        }
      };
      virtual std::string className() const { return "LaDa::MinimizerOp"; }

      bool operator() (EO_OBJECT &_object) 
      {
        object->set_variables( _object.get_variables() );
        minimizer->minimize();
       
        return true;
      }
  }; // class MinimizationOp : public eoMonOp<EO_OBJECT> 

} // endif LaDa

#endif

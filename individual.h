#ifndef _INDIVIDUAL_H_
#define _INDIVIDUAL_H_

#include <eo/eoScalarFitness.h>
#include <opt/opt_function_base.h>
#include <exception>
#include <iostream>
#include <stdexcept>       // std::runtime_error
#include <eo/eoObject.h>      // eoObject
#include <eo/eoPersistent.h>  // eoPersistent
#include <eo/eoOpContainer.h>  

#include<algorithm>
#include<functional>

// functional should be some class derived from 
// opt::Function<TYPE, CONTAINTER>, where CONTAINER supports random
// access iterators
//

namespace LaDa 
{
  template<class FITNESS = eoMinimizingFitness, class FUNCTIONAL = opt::Function<> >
  class Individual : public FUNCTIONAL, public eoObject, public eoPersistent
  {

    public:
      typedef typename FUNCTIONAL::TYPE TYPE;
      typedef typename FUNCTIONAL::CONTAINER CONTAINER;
      typedef typename CONTAINER::iterator CONTAINER_ITERATOR;
      typedef typename CONTAINER::iterator iterator;
      typedef typename CONTAINER::const_iterator CONST_CONTAINER_ITERATOR;
      typedef TYPE AtomType;
      typedef FITNESS Fitness;

    public:
      using FUNCTIONAL :: variables;
      using FUNCTIONAL :: set_variables;
      using FUNCTIONAL :: get_variables;

    private:
      Fitness repFitness;
      int MotU_ref; // ref to master of the universe
      static bool baseline_is_valid;
      bool quantity_is_valid;
      TYPE quantity, baseline;
      unsigned age;
      

    public: 
      Individual() : FUNCTIONAL(), repFitness(Fitness()), 
                     MotU_ref(0), quantity_is_valid (false)
        { variables = new CONTAINER; age = 0; };
      
      Individual(const Individual<FITNESS, FUNCTIONAL> &_indiv ) :
              repFitness(_indiv.repFitness),
              MotU_ref(_indiv.MotU_ref), 
              quantity_is_valid (_indiv.quantity_is_valid),
              quantity(_indiv.quantity), baseline(_indiv.baseline), age(0)
      {
        variables = new CONTAINER(_indiv.variables->size()); 
        *variables = *_indiv.variables;
      };

      ~Individual()
        { if (variables) delete variables;  variables = NULL; };

      void invalidate() { quantity_is_valid = false; }
      void invalidate_baseline()
      { 
        baseline_is_valid = false;
        invalidate();
      } 
      void validate_baseline()
        { baseline_is_valid = true; } 
      bool is_baseline_valid() const
        { return baseline_is_valid; }
      bool invalid() const
        { return not quantity_is_valid; }


      void set_quantity( TYPE _q )
      {
        quantity = _q;
        quantity_is_valid = true;
      }
      TYPE get_quantity() const
        { return quantity; }
      void set_baseline( TYPE _b )
        { baseline = _b; }
      TYPE get_baseline() const
        { return baseline; }

      void  set_age( unsigned _age )
        { age = _age; }
      unsigned  get_age() const
        { return age; }
      void  set_MotU_ref( int _m )
        { MotU_ref = _m; }
      int get_MotU_ref() const
        { return MotU_ref; }

      bool operator<(const Individual<FITNESS, FUNCTIONAL>& _eo2) const
        { return fitness() < _eo2.fitness(); }
      bool operator>(const Individual<FITNESS, FUNCTIONAL>& _eo2) const
        { return !(fitness() <= _eo2.fitness()); }
      void operator=( const Individual<FITNESS, FUNCTIONAL> &_indiv )
      {
        repFitness = _indiv.repFitness;
        MotU_ref   = _indiv.MotU_ref;
        quantity_is_valid = _indiv.quantity_is_valid;
        quantity = _indiv.quantity;
        baseline = _indiv.baseline;
        *variables = *_indiv.variables;
        age = _indiv.age;
      }
      bool operator==( const Individual<FITNESS, FUNCTIONAL> &_indiv )
      {
        unsigned size = variables->size();
        if ( size != _indiv.variables->size() )
          return false;
        if ( size == 0 )
          return true;
        return std::equal( variables->begin(), variables->end(),
                           _indiv.variables->begin() );
      }
        
      Fitness fitness() const
      {
         if (not quantity_is_valid)
           throw std::runtime_error("wtf invalid fitness");
         return repFitness;
      }
      void set_fitness()
      {
        repFitness = quantity - baseline;
        quantity_is_valid = true;
      }

      // eo stuff
      virtual std::string className() const
        { return "Individual<...>"; }
      void printOn(std::ostream &_os) const
      {
        { // EO stuff
          if (invalid()) {
              _os << "INVALID ";
          }
          else
          {
              _os << repFitness << ' ';
          }
        } // EO stuff
        // prints ref
        _os << MotU_ref << " " << variables->size() << " ";
        CONST_CONTAINER_ITERATOR i_var  = variables->begin();
        CONST_CONTAINER_ITERATOR i_last = variables->end();
        for( ; i_var != i_last; ++i_var )
          _os << *i_var << " ";
      }
      void readFrom(std::istream &_is)
      {

        { // EO stuff
          std::string fitness_str;
          int pos = _is.tellg();
          _is >> fitness_str;

          if (fitness_str == "INVALID")
          {
             quantity_is_valid = false;
          }
          else
          {
             quantity_is_valid = true;
             _is.seekg(pos); // rewind
             _is >> repFitness;
          }
        } // EO stuff
        unsigned var_size;
        _is >> var_size;
        variables->resize(var_size);

        CONTAINER_ITERATOR i_var = variables->begin();
        CONTAINER_ITERATOR i_last = variables->end();
        for( ; i_var != i_last; ++i_var )
          _is >> *i_var;
      }

      void print_out( std::ostream &_stream ) const
      {
        CONTAINER_ITERATOR i_var = variables->begin();
        CONTAINER_ITERATOR i_last = variables->end();
        for( ; i_var != i_last; ++i_var )
          _stream << int(*i_var) << " ";
      }

  };


  // redefines operators < and > with respect to FITNESS
  template <class FUNCTIONAL, class FITNESS>
  bool operator<(const Individual<FUNCTIONAL, FITNESS>& _eo1,
                 const Individual<FUNCTIONAL, FITNESS>& _eo2)
  {
      return _eo1.operator<(_eo2);
  }
  template <class FUNCTIONAL, class FITNESS>
  bool operator>(const Individual<FUNCTIONAL, FITNESS>& _eo1,
                 const Individual<FUNCTIONAL, FITNESS>& _eo2)
  {
      return _eo2.operator<(_eo1);
  }

  template<class FITNESS, class FUNCTIONAL> 
  bool Individual<FITNESS, FUNCTIONAL> :: baseline_is_valid = true; 
} // endif LaDa

#endif

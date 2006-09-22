#ifndef _INDIVIDUAL_H_
#define _INDIVIDUAL_H_

#include <eo/EO.h>
#include <eo/eoScalarFitness.h>
#include <opt/opt_function_base.h>
#include <exception>
#include <iostream>

// functional should be some class derived from 
// opt::Function<TYPE, CONTAINTER>, where CONTAINER supports random
// access iterators
//

namespace LaDa 
{
  template<class FUNCTIONAL = opt::Function<double>, class FITNESS = eoMaximizingFitness> 
  class Individual : public FUNCTIONAL,
                     public EO<FITNESS> 
  {

    public:
      typedef typename FUNCTIONAL::TYPE TYPE;
      typedef typename FUNCTIONAL::CONTAINER CONTAINER;
      typedef typename FUNCTIONAL::CONTAINER_ITERATOR CONTAINER_ITERATOR;
      typedef typename FUNCTIONAL::CONST_CONTAINER_ITERATOR CONST_CONTAINER_ITERATOR;
      typedef typename FUNCTIONAL::TYPE AtomType;
      using EO<FITNESS>::invalid;
      using EO<FITNESS>::invalidate;

    protected:
      using FUNCTIONAL::variables;

    private:
      int MotU_ref; // ref to master of the universe
      static bool baseline_is_valid;
      bool quantity_is_valid;
      TYPE quantity, baseline;
      CONTAINER chromosome;

    public: 
      Individual( int _i = -1 )
        { variables = &chromosome; MotU_ref = _i; }

      void invalidate_quantity()
      {
        quantity_is_valid = false;
        invalidate();
      } 
      void invalidate_baseline()
      { 
        baseline_is_valid = false;
        invalidate();
      } 
      void validate_baseline()
        { baseline_is_valid = true; } 
      void validate_quantity()
        { quantity_is_valid = true; } 
      bool is_baseline_valid() const
        { return baseline_is_valid; }
      bool is_quantity_valid() const
        { return quantity_is_valid; }
      bool invalidate_quantity() const
      {
        quantity_is_valid = false; 
        invalidate();
      }


      void set_quantity( TYPE _q )
      {
        quantity = _q;
        quantity_is_valid = true;
        invalidate();
      }
      void set_baseline( TYPE _b )
        { baseline = _b; }

      void  set_MotU_ref( int _m )
        { MotU_ref = _m; }
      int get_MotU_ref() const
        { return MotU_ref; }

      // eo stuff
      virtual std::string className() const
        { return "Individual<...>"; }
      void printOn(std::ostream &_os) const
      {
        if (not variables or MotU_ref < 0) 
          return;

        // prints ref
        _os << MotU_ref << " " << variables << " ";
        CONST_CONTAINER_ITERATOR i_var  = variables->begin();
        CONST_CONTAINER_ITERATOR i_last = variables->end();
        for( ; i_var != i_last; ++i_var )
          _os << *i_var << " ";
      }

      void set_fitness()
        { fitness( quantity - baseline ); }

      void readFrom(std::istream &_is)
      {
        EO<FITNESS>::readFrom(_is);

        unsigned var_size;
        _is >> var_size;
        if( not variables )
          variables = new CONTAINER;
        variables->resize(var_size);

        CONTAINER_ITERATOR i_var = variables->begin();
        CONTAINER_ITERATOR i_last = variables->end();
      }
      TYPE& operator[](size_t n)
      {
        if ( n > variables->size() )
          throw std::out_of_range("out of range in Individual<...> :: TYPE& operator[](size_t n) ");
            
        return *(variables->begin() + n ); 
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

  template<class FUNCTIONAL, class FITNESS> 
  bool Individual<FUNCTIONAL, FITNESS> :: baseline_is_valid = true; 
} // endif LaDa

#endif

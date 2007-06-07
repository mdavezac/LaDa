// Implements a base template class for functions minimizable by
// opt_minimize. It contains the following members
//
//   CONTAINER *variables: container with the function's variable
//   typename t_Container :: iterator iterator, needed
//     for iterating over variables
//   typename t_Container :: const_iterator const_iterator,
//     same but constant  
//   Base, constructor, doesn't do anything;
//   ~Base, destructor, doesn't do anything;
//   VARIABLE_TYPE_ITERATOR begin() const, returns iterator to
//       beginning of container
//   VARIABLE_TYPE_ITERATOR end() const, returns iterator to
//       end of container
//   void reserve(const types::t_unsigned _nb), sets size of container
//   types::t_unsigned size() const, returns size of container
//    _T get_concentration() const, average over container
//
// You should provide the following member functions
// _T evaluate() const; evaluates the function
// void evaluate_gradient(_T *gradient) const; evaluates the gradient
//    of  the function and stores it in gradient
// _T evaluate_gradient(_T *gradient) const; evaluates funtion and its gradient
// 
#ifndef _OPT_FUNCTION_BASE_H_
#define _OPT_FUNCTION_BASE_H_

#include <iostream>
#include <algorithm>
#include <numeric>
#include <vector>
#include <stdexcept>
#include <opt/types.h>


namespace function
{

  template<class _TYPE=types::t_real, class _CONTAINER = std::vector<_TYPE> >
  class Base
  {
    public:
      typedef _TYPE t_Type;
      typedef _CONTAINER t_Container;
      typedef typename t_Container :: iterator iterator;
      typedef typename t_Container :: const_iterator const_iterator;

    protected:
      t_Container *variables;
      bool does_own_variables;

    public:
      Base() : variables(NULL), does_own_variables(false) {}
      Base( types::t_int nb ) : variables(NULL), does_own_variables(false) 
        { resize(nb); }
      Base   ( const Base<t_Type, t_Container> &_f )
           : variables( _f.variables ), does_own_variables( false ) {}
      virtual ~Base() { if ( variables and does_own_variables ) delete variables; }
      
      iterator begin() 
        { return variables->begin(); }
      iterator end() 
        { return variables->end(); }
      const_iterator begin() const
        { return variables->begin(); }
      const_iterator end() const 
        { return variables->end(); }
      void reserve(const types::t_unsigned _nb) 
      {
        if ( not variables )
        {
          variables = new t_Container;
          does_own_variables = true; 
        }
        variables->reserve(_nb); 
      }
      virtual void resize(const types::t_unsigned _nb) 
      {
        if ( not variables )
        {
          variables = new t_Container;
          does_own_variables = true; 
        }
        variables->resize(_nb, t_Type(0) ); 
      }
      types::t_unsigned size() const 
      { return variables->size(); }
      t_Type get_concentration() const
      {
         return ( std::accumulate(variables->begin(), variables->end(), 0.0) 
                  / ( (t_Type) variables->size() ) );
      }
      virtual void set_variables( t_Container* const var ) 
       { variables = var; does_own_variables = false; }
      virtual t_Container* get_variables() const
       { return variables; }
      void get_variables( t_Container *var ) const
      {
         if ( var->size() < variables->size() )
           var->resize( variables->size(), 1.0 );
         std::copy( variables->begin(), variables->end(), var->begin() );
      }
      void copy_variables( t_Container *var ) const
      {
         if ( variables->size() < var->size() )
         { 
           std::cerr << " error while copying variables in function::Base::copy_variables "
                     << std::endl;
           exit(0);
         }
         std::copy( var->begin(), var->end(), variables->begin() );
      }

      void destroy_variables()
      {
        if ( does_own_variables and variables )
          delete variables;
        variables = NULL;
        does_own_variables = false;
      }

      t_Type& operator[](size_t n)
      {
        if ( n > variables->size() )
          throw std::out_of_range("out of range in function::Base<...> :: t_Type& operator[](size_t n) ");
            
        return *(variables->begin() + n ); 
      }

      virtual bool init() {return true;}; // sometimes, we should pass stuff from somewhere to variables
        
      virtual t_Type evaluate() = 0; 
      virtual void evaluate_gradient( t_Type* const _i_grad ) = 0;
      virtual t_Type evaluate_one_gradient( types::t_unsigned _pos) = 0;
      virtual t_Type evaluate_with_gradient( t_Type* const _i_grad ) = 0;

      virtual bool is_taboo() const { return false; } // for minimizers with taboos
  };

}

#endif // _OPT_FUNCTION_BASE_H_

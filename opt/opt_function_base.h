//
//  Version: $Id$
//
#ifndef _OPT_FUNCTION_BASE_H_
#define _OPT_FUNCTION_BASE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <algorithm>
#include <numeric>
#include <vector>
#include <stdexcept>

#include "opt/types.h"


//! \brief Brings together an abstract functional base class and functors for that class
//! \details This base class is expected to be derived from by true functionals
//! such as Vff::Functional. It declares an interface with which minimizer
//! objects can use to do their own job.
namespace function
{

  //! \brief Abstract base class for all functionals
  //! \details The object of this class is to declare an interface with which
  //! objects in namespace minimizer can work.
  //! Base is templated as follows
  //! \param _TYPE is the type of the variable. If Base is \f$f(x)\f$, then  \a
  //!        _TYPE is the type of \f$x\f$ (with \f$x\f$ a scalar).
  //! \param _CONTAINER is a container, say a vector which is \f$x\f$
  //! 
  //! The variables of the functional are declared in Base a pointer to a
  //! container, Base::variables. It is the values within this container which
  //! a minimizer will change to get to the minimum. Base::variables may or may
  //! not be owned by an instance of Base. This class takes care of that, and
  //! destroys Base::variables only the particular instance of Base owns the
  //! pointer to Base::variables.
  //! For summing functionals and such operations, see Function::Plus, Function::Minus
  //
  //! \par Abstract functions
  //! A number of abstract functions are declared which all deal with evaluating the functional
  //! - Base::evaluate() returns the value of the functional 
  //! - Base::evaluate_with_gradient() returns both the value and the gradients of the functional 
  //! - Base::evaluate_gradient() returns the gradients of the functional 
  //! - Base::evaluate_one_gradient() returns the gradient in a single direction
  //! .
  //! \todo use Traits::Function instead of the rather pathetic _TYPE and _CONTAINER
  template<class _TYPE=types::t_real, class _CONTAINER = std::vector<_TYPE> >
  class Base
  {
    public:
      typedef _TYPE t_Type; //!< Type of the (scalar) variable
      typedef _CONTAINER t_Container; //!< Type of the full vectorial variable
      typedef typename t_Container :: iterator iterator; //!< Iterator type
      typedef typename t_Container :: const_iterator const_iterator;//!< Constant Iteratot type

    protected:
      t_Container *variables; //!< Pointer to the variables of the functional
      bool does_own_variables;//!< Whether this instance of Base owns Base::variables

    public:
      //! \brief Constructor and Initializer
      //! \details can assign an external container as Base::variables
      Base( t_Container *_vars = NULL ) : variables(_vars), does_own_variables(false) {}
      //! \brief Constructor and Initializer
      //! \details creates a container Base::variables. After call,
      //! Base::variables should be a newly created container owned by this
      //! instance of Base.
      Base( types::t_int nb ) : variables(NULL), does_own_variables(false) 
        { resize(nb); }
      //! Copy Constructor
      Base   ( const Base<t_Type, t_Container> &_f )
           : variables( _f.variables ), does_own_variables( false ) {}
      //! \brief Deconstructor
      //! \details destroys Base::variables if an only if Base::variables is
      //! non-NULL and owned by this particular instance of Base, according to
      //! Base::does_own_variables.
      virtual ~Base() { if ( variables and does_own_variables ) delete variables; }
      
      //! \brief Returns an iterator to the start of variables. 
      //! \details Does not check if variables is non-zero
      iterator begin() 
        { return variables->begin(); }
      //! \brief Returns an iterator to the end of variables. 
      //! \details Does not check if variables is non-zero
      iterator end() 
        { return variables->end(); }
      //! \brief Returns a constant iterator to the beginning of variables. 
      //! \details Does not check if variables is non-zero
      const_iterator begin() const
        { return variables->begin(); }
      //! \brief Returns a constant iterator to the end of variables. 
      //! \details Does not check if variables is non-zero
      const_iterator end() const 
        { return variables->end(); }
      //! \brief applies reserve to Base::variables
      //! \details If Base::variables is NULL, then it is first created. If it
      //! is created, then it is owned by this particular instance of Base.
      void reserve(const types::t_unsigned _nb);
      //! \brief applies resize to Base::variables
      //! \details If Base::variables is NULL, then it is first created. If it
      //! is created, then it is owned by this particular instance of Base.
      virtual void resize(const types::t_unsigned _nb);
      //! \brief returns the size of Base::variables
      //! \details Does not check if variables is non-zero
      types::t_unsigned size() const
        { return variables->size(); }
      //! \brief returns the mean value of Base::variables
      //! \details Does not check if variables is non-zero
      t_Type get_concentration() const;
      //! \brief Sets Base::variables to \a _var
      //! \details Does not destory variables if it is non-NULL before
      //! assignement. Does not check tha \a _var is valid. After this call,
      //! Base::variables is not owned by this particular instance of Base.
      virtual void set_variables( t_Container* const _var ) 
       { variables = _var; does_own_variables = false; }
      //! \brief Returns the value of Base::variables
      virtual t_Container* get_variables() const
       { return variables; }
      //! \brief Copies Base::variables into \a _var
      void get_variables( t_Container *_var ) const;
      //! \brief Copies \a _var into Base::variables 
      void copy_variables( t_Container *_var ) const;

      //! \brief Destroys Base::variables if is both valid and owned. 
      //! \details After call, Base::variables is invalid.
      void destroy_variables();

      //! \brief Returns component \a _n of Base::variables
      t_Type& operator[](size_t n);

      //! \brief Should container whatever initialization may be necessary before minimization
      virtual bool init() {return true;}; 
        
      //! \brief should return the value of the functional for the current Base::variables
      virtual t_Type evaluate() = 0; 
      //! \brief should return the gradient of the functional for the current Base::variables
      virtual void evaluate_gradient( t_Type* const _i_grad ) = 0;
      //! \brief should return the gradient of the functional in direction \a _pos 
      //! for the current Base::variables
      virtual t_Type evaluate_one_gradient( types::t_unsigned _pos) = 0;
      //! \brief should return the value and the gradient of the functional
      //! for the current Base::variables
      virtual t_Type evaluate_with_gradient( t_Type* const _i_grad ) = 0;

      //! \brief Should return true if the current state of Base::variables is Taboo
      //! \sa minimizers::VA
      virtual bool is_taboo() const { return false; } // for minimizers with taboos

      //! \brief For expensive functionals, this member function helps to make sure the
      //! functional is computed only once
      virtual void invalidate() {}
      //! \brief For expensive functionals, this member function helps to make sure the
      //! functional is computed only once
      virtual bool is_valid() { return true; }
  };


  template<class _TYPE, class _CONTAINER >
    inline void Base<_TYPE,_CONTAINER> :: reserve(const types::t_unsigned _nb) 
    {
      if ( not variables )
      {
        variables = new t_Container;
        does_own_variables = true; 
      }
      variables->reserve(_nb); 
    }
  template<class _TYPE, class _CONTAINER >
    inline void Base<_TYPE,_CONTAINER> :: resize(const types::t_unsigned _nb) 
    {
      if ( not variables )
      {
        variables = new t_Container;
        does_own_variables = true; 
      }
      variables->resize(_nb, t_Type(0) ); 
    }
  template<class _TYPE, class _CONTAINER >
    inline Base<_TYPE,_CONTAINER> :: t_Type Base<_TYPE,_CONTAINER> :: get_concentration() const
    {
       return ( std::accumulate(variables->begin(), variables->end(), 0.0) 
                / ( (t_Type) variables->size() ) );
    }
  template<class _TYPE, class _CONTAINER >
    inline void Base<_TYPE,_CONTAINER> :: get_variables( t_Container *_var ) const
    {
       if ( _var->size() < variables->size() )
         _var->resize( variables->size(), 1.0 );
       std::copy( variables->begin(), variables->end(), _var->begin() );
    }
  template<class _TYPE, class _CONTAINER >
    inline void Base<_TYPE,_CONTAINER> :: copy_variables( t_Container *var ) const
    {
       if ( variables->size() < _var->size() )
       { 
         std::cerr << " error while copying variables in function::Base::copy_variables "
                   << std::endl;
         exit(0);
       }
       std::copy( _var->begin(), _var->end(), variables->begin() );
    }

  template<class _TYPE, class _CONTAINER >
    inline void Base<_TYPE,_CONTAINER> ::  destroy_variables()
    {
      if ( does_own_variables and variables )
        delete variables;
      variables = NULL;
      does_own_variables = false;
    }

  template<class _TYPE, class _CONTAINER >
    inline Base<_TYPE,_CONTAINER>::t_Type& Base<_TYPE,_CONTAINER> :: operator[](size_t n)
    {
      if ( n > variables->size() )
        throw std::out_of_range("out of range in function::Base<...> :: t_Type& operator[](size_t n) ");
          
      return *(variables->begin() + n ); 
    }

}

#endif // _OPT_FUNCTION_BASE_H_

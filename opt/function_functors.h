//
//  Version: $Id$
//
#ifndef _OPT_FUNCTORS_H_
#define _OPT_FUNCTORS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<iostream>

#include "function_base.h"
#include "opt/types.h"

namespace function
{
  //! \brief Implements a wrapper around a function::Base to do such things a
  //!        mutlipling that function::Base by itself.
  //! \details This is mostly a general base class, used by, for instance,
  //!          Quadratic below.
  template<typename T_FUNCTIONAL>
  class Unary_Functor : public Base<typename T_FUNCTIONAL::t_Type,
                                    typename T_FUNCTIONAL::t_Container>
  {

    public:
      //! The type of the functional
      typedef T_FUNCTIONAL t_Functional;
      //! see function::Base::t_Type
      typedef typename t_Functional :: t_Type t_Type;
      //! see function::Base::t_Container
      typedef typename t_Functional :: t_Container t_Container;

    protected:
      //! Type of the base classa
      typedef Base<t_Type, t_Container> t_Base;

    protected:
      using t_Base :: variables;
      using t_Base :: does_own_variables;

    public:
      using Base<t_Type, t_Container> :: size;

    protected:
      //! A pointer to the functional 
      t_Functional *functional;

    public:     
      //! Constructor and Initializer. 
      explicit 
        Unary_Functor   (const t_Functional* _func)
                      : t_Base(_func), functional(_func) {}
      //! Copy Constructor. 
      Unary_Functor   (const Unary_Functor< T_FUNCTIONAL > &_f)
                    : t_Base(_f), functional(_f) {}
      //! Constructor.
      Unary_Functor() : t_Base(), functional(NULL) {}
                   
      //! Destructor.
      virtual ~Unary_Functor() {};
      
      //! Sets the functional, and optionally the variables.
      void set_functional( t_Functional *_func, bool set_variables_to_func = true );
      //! Returns the functional
      t_Functional* get_functional() const { return functional; };
    
    public:
      //! Sets the variables
      virtual void set_variables( t_Container* const var );
      //! Destroys the variables
      virtual void destroy_variables();
      //! Returns a pointer to the variables
      virtual t_Container* get_variables() const { return  variables; }

      //! Resizes the number of variables
      virtual void resize(const types::t_unsigned _nb);

      //! Calls the init() member of the functional
      bool init() { return functional->init(); }

      //! Calls the is_taboo() member of the functional
      bool is_taboo() const { return functional->is_taboo(); }
  };
 
  //! \brief Implements a wrapper around two function::Base to do such things a
  //!        mutlipling these function::Base.
  //! \details This is mostly a general base class, used by, for instance,
  //!          Plus or Minus below. There is no assurance that binary functors
  //!          here or below will work if both functionals take different
  //!          %types of arguments and/or return different %types of values.
  //!          The variables used are those of the first functional. The second
  //!          functional is made to point to the variables of the first functional.
  template<typename T_FUNC1, typename T_FUNC2>
  class Binary_Functor : public Base<typename T_FUNC1::t_Type, typename T_FUNC2::t_Container>
  {
    public:
      //! The type of the first functional
      typedef T_FUNC1 t_Functional1;
      //! The type of the second functional
      typedef T_FUNC2 t_Functional2;
      //! see function::Base::t_Type
      typedef typename t_Functional1 :: t_Type t_Type;
      //! see function::Base::t_Container
      typedef typename t_Functional1 :: t_Container t_Container;

    protected:
      //! Type of the base classa
      typedef Base<t_Type, t_Container> t_Base;

    protected:
      using Base<t_Type, t_Container> :: variables;
      using Base<t_Type, t_Container> :: does_own_variables;

    public:
      using Base<t_Type, t_Container> :: size;

    protected:
      //! A pointer to the first functional 
      t_Functional1 *functional1;
      //! A pointer to the second functional 
      t_Functional2 *functional2;

    public:  
      //! Constructor and Initializer.
      explicit Binary_Functor   ( const t_Functional1 *_functional1,
                                  const t_Functional2 *_functional2 );
      //! Copy Constructor.
      Binary_Functor   (const Binary_Functor<t_Functional1, t_Functional2> &_f) 
                     : t_Base( _f ), functional1( _f.functional1 ),
                       functional2( _f.functional2 ) {}
      //! Constructor
      Binary_Functor() : t_Base(), functional1(NULL), functional2(NULL) {}
                   
      //! Destructor
      virtual ~Binary_Functor() {};
      
      //! Sets the first functional.
      void set_functional1( t_Functional1 * const _functional1,
                            bool set_variables_to_obj = true );
      //! Sets the second functional.
      void set_functional2( t_Functional2 *const _functional2 );
      //! Returns  a pointer to the first functional
      t_Functional1* get_functional1() const { return functional1; };
      //! Returns  a pointer to the second functional
      t_Functional2* get_functional2() const { return functional2; };
    
    public:
      //! Sets the variables
      virtual void set_variables( t_Container* const var );
      //! Destroys the variables
      virtual void destroy_variables();
      //! Returns a pointer to the variables
      virtual t_Container* get_variables() const { return  variables; }

      //! Resizes the number of variables
      virtual void resize(const types::t_unsigned _nb);

      //! Calls the init() member of both functionals
      virtual bool init()
       { return functional1->init() and functional2->init(); }

      //! Calls the is_taboo() member of both functional
      virtual bool is_taboo() const
       { return functional1->is_taboo() or functional2->is_taboo(); }
  };

  //! Does the difference of two functionals
  template<typename t_Functional1, typename t_Functional2>
  class Minus : public Binary_Functor<t_Functional1, t_Functional2>
  {
    public:
      //! see function::Base::t_Type
      typedef typename Binary_Functor<t_Functional1, t_Functional2> :: t_Type t_Type;
      //! see function::Base::t_Container
      typedef typename Binary_Functor<t_Functional1, t_Functional2> :: t_Container t_Container;

    protected:
      //! Type of the base class
      typedef Binary_Functor< t_Functional1, t_Functional2> t_Base;

    protected:
      using t_Base :: functional1;
      using t_Base :: functional2;
      using t_Base :: variables;
      using t_Base :: does_own_variables;
    public:
      using t_Base :: size;


    public:     
      //! Constructor and Initializer
      explicit Minus  (const t_Functional1 *_f1, const t_Functional2 *_f2)
                     : t_Base( _f1, _f2 ) {}
      //! Constructor
      Minus() : t_Base() {};
      //! Destructor
      virtual ~Minus() {};
      
    public:
      //! Returns the difference of the two functionals
      virtual t_Type evaluate() 
        { return functional1->evaluate() - functional2->evaluate(); }
      //! Computes the gradient of the difference of the two functionals
      virtual void evaluate_gradient( t_Type* const _i_grad);
      //! Computes the gradient of the difference of the two functionals in direction \a _pos
      virtual t_Type evaluate_one_gradient( types::t_unsigned _pos);
      //! Computes the difference of the two functionals and its gradient
      virtual t_Type evaluate_with_gradient( t_Type* const _i_grad);
  };

  //! Does the sum of two functionals
  template<typename t_Functional1, typename t_Functional2>
  class Plus : public Binary_Functor<t_Functional1, t_Functional2>
  {
    public:
      //! see function::Base::t_Type
      typedef typename Binary_Functor<t_Functional1, t_Functional2> :: t_Type t_Type;
      //! see function::Base::t_Container
      typedef typename Binary_Functor<t_Functional1, t_Functional2> :: t_Container t_Container;

    protected:
      //! Type of the base class
      typedef Binary_Functor< t_Functional1, t_Functional2> t_Base;

    protected:
      using t_Base :: functional1;
      using t_Base :: functional2;
      using t_Base :: variables;
      using t_Base :: does_own_variables;
    public:
      using t_Base :: size;

    public:   
      //! Constructor and Initializer
      explicit Plus  (const t_Functional1 *_f1, const t_Functional2 *_f2)
                     : t_Base( _f1, _f2 ) {}
      //! Constructor
      Plus() : t_Base() {};
      //! Destructor
      virtual ~Plus(){};
      
      // required minimizer:: behavior
    public:
      //! Returns the sum of the two functionals
      virtual t_Type evaluate() 
        { return (functional1->evaluate() + functional2->evaluate()); }
      //! Computes the gradient of the sum of the two functionals
      virtual void evaluate_gradient( t_Type* const _i_grad);
      virtual t_Type evaluate_one_gradient( types::t_unsigned _pos);
      //! Computes the sum of the two functionals and its gradient
      virtual t_Type evaluate_with_gradient( t_Type* const _i_grad );
  };

  //! Does the square of a functional, \f$ (F(x) - y )^2 \f$.
  template<typename T_FUNCTIONAL>
  class Quadratic : public Unary_Functor< T_FUNCTIONAL >
  {
    public:
      //! Type of the functional
      typedef T_FUNCTIONAL t_Functional;
      //! see function::Base::t_Type
      typedef typename t_Functional :: t_Type t_Type;
      //! see function::Base::t_Container
      typedef typename t_Functional :: t_Container t_Container;

    protected:
      //! Type of the base class
      typedef Unary_Functor< t_Functional > t_Base;

    protected:
      using t_Base :: functional;
      using t_Base :: size;
      t_Type center; //!< reference shift

    public:
      //! Constructor and initializer.
      explicit Quadratic   (const t_Functional *_func, const t_Type &_center = t_Type(0) )
                         : t_Base(_func), center(_center) {}
      //! Constructor and initializer.
      Quadratic   ( const t_Type &_center = t_Type(0) )
                : t_Base(), center(_center) {};
      //! Destructor
      virtual ~Quadratic(){};
      
    public:
      //! Returns \f$ (F(x) - y )^2 \f$.
      virtual t_Type evaluate() 
        { t_Type value = functional->evaluate() - center; return value * value; }
      //! Returns \f$ (F(x) - y )^2 \f$ and its gradient with respect to \e x.
      virtual t_Type evaluate_with_gradient( t_Type* const _i_grad);
      //! Returns the gradient of \f$ (F(x) - y )^2 \f$ with respect to \e x.
      virtual void evaluate_gradient( t_Type* const _i_grad ) 
        { evaluate_with_gradient(_i_grad); }
      //! Returns the gradient of \f$ (F(x) - y )^2 \f$ with respect to \e x[\a _pos].
      virtual t_Type evaluate_one_gradient( types::t_unsigned _pos);
  };



  template<typename T_FUNCTIONAL>
    inline void Unary_Functor<T_FUNCTIONAL> :: set_functional( t_Functional *_func, 
                                                               bool set_variables_to_func )
    { 
      functional = _func; 
      if ( set_variables_to_func )
      {
        variables = _func->get_variables();
        does_own_variables = false; 
      }
    }
  template<typename T_FUNCTIONAL>
    inline void Unary_Functor<T_FUNCTIONAL> :: set_variables( t_Container* const var ) 
    {
      variables = var;
      does_own_variables = false;
      functional->set_variables(var);
    }
  template<typename T_FUNCTIONAL>
    inline void Unary_Functor<T_FUNCTIONAL> :: destroy_variables()
    {
      functional->destroy_variables();
      variables = NULL;
      does_own_variables = false;
    }
  template<typename T_FUNCTIONAL>
    inline void Unary_Functor<T_FUNCTIONAL> :: resize(const types::t_unsigned _nb) 
    {
      functional->resize(_nb);
      variables = functional->get_variables();
      does_own_variables = false;
    }




  template<typename T_FUNC1, typename T_FUNC2>
    Binary_Functor<T_FUNC1, T_FUNC2> :: Binary_Functor   ( const t_Functional1 *_functional1,
                                                           const t_Functional2 *_functional2 )
                                                       : functional1(_functional1),
                                                         functional2(_functional2)
    { 
      variables = functional1->get_variables();
      does_own_variables = false;
      functional2->set_variables(variables); 
    }
    
  template<typename T_FUNC1, typename T_FUNC2>
    inline void Binary_Functor<T_FUNC1, T_FUNC2> :: 
      set_functional1( t_Functional1 * const _functional1,
                       bool set_variables_to_obj )
      { 
        functional1 = _functional1; 
        if ( set_variables_to_obj )
        {
          variables = functional1->get_variables();
          does_own_variables = false;
        }
        else
          functional1->set_variables( variables );
      
        if ( set_variables_to_obj and functional2 )
          functional2->set_variables( variables );
      }
  template<typename T_FUNC1, typename T_FUNC2>
    inline void Binary_Functor<T_FUNC1, T_FUNC2> ::  
      set_functional2( t_Functional2 *const _functional2 )
      { 
        functional2 = _functional2; 
        functional2->set_variables(variables);
      }
  template<typename T_FUNC1, typename T_FUNC2>
    inline void Binary_Functor<T_FUNC1, T_FUNC2> :: set_variables( t_Container* const var ) 
    {
      variables = var;
      does_own_variables = false;
      functional1->set_variables(var);
      functional2->set_variables(var);
    }
  template<typename T_FUNC1, typename T_FUNC2>
    inline void Binary_Functor<T_FUNC1, T_FUNC2> :: destroy_variables()
    {
      if (does_own_variables and variables) delete variables;
      variables = NULL; does_own_variables = false;
      functional1->set_variables(NULL);
      functional2->set_variables(NULL);
    }
  template<typename T_FUNC1, typename T_FUNC2>
    inline void Binary_Functor<T_FUNC1, T_FUNC2> :: resize(const types::t_unsigned _nb) 
    {
      functional1->resize(_nb);
      variables = functional1->get_variables();
      does_own_variables = false;
      functional2->set_variables(variables);
    }




  template<typename t_Functional1, typename t_Functional2>
    inline void Minus<t_Functional1, t_Functional2> ::
      evaluate_gradient( t_Type* const _i_grad) 
      {
        types::t_unsigned N = size();
        t_Type *grad = new t_Type[N];
        if ( !grad )
        { 
          std::cerr << "Could not allocate memory in function::Minus"
                    << std::endl;
          exit(1);
        }
        functional1->evaluate_gradient(_i_grad);
        functional2->evaluate_gradient(grad);
      
        t_Type *ptr_grad1, *ptr_grad2, *ptr_grad_last;
        ptr_grad1 = grad;
        ptr_grad2 = _i_grad;
        ptr_grad_last = _i_grad + N;
        for ( ; ptr_grad2 != ptr_grad_last; ++ptr_grad1, ++ptr_grad2)
          *ptr_grad2 -= *ptr_grad1;
        
        delete[] grad;
      }
  template<typename t_Functional1, typename t_Functional2>
    inline typename Minus<t_Functional1, t_Functional2> :: t_Type
      Minus<t_Functional1, t_Functional2> :: evaluate_one_gradient( types::t_unsigned _pos) 
      {
        return   functional1->evaluate_one_gradient(_pos)
               + functional2->evaluate_one_gradient(_pos);
      }
  template<typename t_Functional1, typename t_Functional2>
    inline typename Minus<t_Functional1, t_Functional2> :: t_Type
      Minus<t_Functional1, t_Functional2> :: evaluate_with_gradient( t_Type* const _i_grad)
      {
        types::t_unsigned N = size();
        t_Type *grad = new t_Type[N];
        if ( !grad )
        { 
          std::cerr << "Could not allocate memory in function::Minus"
                    << std::endl;
          exit(1);
        }
        t_Type value  = functional1->evaluate_with_gradient(_i_grad);
        value -= functional2->evaluate_with_gradient(grad);

        t_Type *ptr_grad1, *ptr_grad2, *ptr_grad_last;
        ptr_grad1 = grad;
        ptr_grad2 = _i_grad;
        ptr_grad_last = _i_grad + N;
        for ( ; ptr_grad2 != ptr_grad_last; ++ptr_grad1, ++ptr_grad2)
          *ptr_grad2 -= *ptr_grad1;
        
        delete[] grad;

        return value;
      }




  template<typename t_Functional1, typename t_Functional2>
    inline void Plus<t_Functional1, t_Functional2> ::
      evaluate_gradient( t_Type* const _i_grad) 
      {
        types::t_unsigned N = size();
        t_Type *grad = new t_Type[N];
        if ( !grad )
        { 
          std::cerr << "Could not allocate memory in function::Plus"
                    << std::endl;
          exit(1);
        }
        functional1->evaluate_gradient(_i_grad);
        functional2->evaluate_gradient(grad);

        t_Type *ptr_grad1, *ptr_grad2, *ptr_grad_last;
        ptr_grad1 = grad;
        ptr_grad2 = _i_grad;
        ptr_grad_last = _i_grad + N;
        for ( ; ptr_grad2 != ptr_grad_last; ++ptr_grad1, ++ptr_grad2)
          *ptr_grad2 += *ptr_grad1;
        
        delete[] grad;
      }
  template<typename t_Functional1, typename t_Functional2>
    inline typename Plus<t_Functional1, t_Functional2> :: t_Type
      Plus<t_Functional1, t_Functional2> :: evaluate_one_gradient( types::t_unsigned _pos) 
      {
        return   functional1->evaluate_one_gradient(_pos)
               + functional2->evaluate_one_gradient(_pos);
      }
  template<typename t_Functional1, typename t_Functional2>
    inline typename Plus<t_Functional1, t_Functional2> :: t_Type
      Plus<t_Functional1, t_Functional2> :: evaluate_with_gradient( t_Type* const _i_grad)
      {
        types::t_unsigned N = size();
        t_Type *grad = new t_Type[N];
        if ( !grad )
        { 
          std::cerr << "Could not allocate memory in function::Plus"
                    << std::endl;
          exit(1);
        }
        t_Type value  = functional1->evaluate_with_gradient(_i_grad);
        value += functional2->evaluate_with_gradient(grad);

        t_Type *ptr_grad1, *ptr_grad2, *ptr_grad_last;
        ptr_grad1 = grad;
        ptr_grad2 = _i_grad;
        ptr_grad_last = _i_grad + N;
        for ( ; ptr_grad2 != ptr_grad_last; ++ptr_grad1, ++ptr_grad2)
          *ptr_grad2 += *ptr_grad1;
        
        delete[] grad;

        return value;
      }




  template<typename T_FUNCTIONAL>
    inline typename Quadratic<T_FUNCTIONAL> :: t_Type 
      Quadratic<T_FUNCTIONAL> :: evaluate_with_gradient( t_Type* const _i_grad) 
      {
        t_Type value = functional->evaluate_with_gradient(_i_grad) - center;
      
        t_Type grad_factor = value * 2.0;
        t_Type *ptr_grad, *ptr_grad_last;
        ptr_grad = _i_grad;
        ptr_grad_last = _i_grad + size();
        for ( ; ptr_grad != ptr_grad_last; ++ptr_grad)
          *ptr_grad *= (*ptr_grad) * grad_factor;
        
        return value * value;
      }
  template<typename T_FUNCTIONAL>
    inline typename Quadratic<T_FUNCTIONAL> :: t_Type 
      Quadratic<T_FUNCTIONAL> :: evaluate_one_gradient( types::t_unsigned _pos) 
      {
        t_Type value = functional->evaluate() - center;
        return  2.0 * functional->evaluate_one_gradient(_pos) * value;
      }
} // namespace function
#endif // _OPT_FUNCTORS_H

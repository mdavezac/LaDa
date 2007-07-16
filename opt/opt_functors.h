// Implements binary functors for operations OP="-" and OP="+"
// for use with minimizer::?? interface
// both classes should have the same variable arrays (ie in mem)
// this works for me...
//
// the following member functions are implemented:
//   double evaluate():
//       returns t_Functional1 OP t_Functional2 with current variables
//   void evaluate_gradient(double* ):
//       stores the gradient of t_Functional1 OP t_Functional2 with current variables
//       into double* input array
//   double evaluate_with_gradient(double* ):
//       stores the gradient of t_Functional1 OP t_Functional2 with current variables
//       into double* input array
//       and returns the t_Functional1 OP t_Functional2 
//   void set_t_Functional1( t_Functional1* ):
//       sets pointer to t_Functional1 object;
//   void set_t_Functional2( t_Functional2* ):
//       sets pointer to t_Functional2 object;
//
//
//
//
// No debugging capability is offered by the functors
// consistency (number of variables, pointers...) 
// should be checked by you!


#ifndef _OPT_FUNCTORS_H_
#define _OPT_FUNCTORS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<iostream>

#include "opt_function_base.h"
#include "opt/types.h"

namespace function
{
  template<typename T_FUNCTIONAL>
  class Unary_Functor : public Base<typename T_FUNCTIONAL::t_Type,
                                    typename T_FUNCTIONAL::t_Container>
  {
    public:
      typedef T_FUNCTIONAL t_Functional;
      typedef typename t_Functional :: t_Type t_Type;
      typedef typename t_Functional :: t_Container t_Container;
      using Base<t_Type, t_Container> :: size;

    protected:
      using Base<t_Type, t_Container> :: variables;
      using Base<t_Type, t_Container> :: does_own_variables;

    protected:
      t_Functional *functional;

    public:     // t_Functional1 and t_Functional2 should exist, so doesn't call
      explicit  // constructors and destructors
        Unary_Functor(const t_Functional* _func) : functional(_func)
          { variables = _func->get_variables(); does_own_variables = false; }
      Unary_Functor  (const Unary_Functor<t_Functional> &_f) 
                    : Unary_Functor( _f.func ) {};
      Unary_Functor() : functional(NULL)
        { variables = NULL; does_own_variables = false; }
                   
      virtual ~Unary_Functor() {};
      
      void set_functional( t_Functional *_func, bool set_variables_to_func = true )
      { 
        functional = _func; 
        if ( set_variables_to_func )
        {
          variables = _func->get_variables();
          does_own_variables = false; 
        }
      }
      t_Functional* get_functional() const { return functional; };
    
      // required minimizer::?? behavior
    public:
      virtual void set_variables( t_Container* const var ) 
      {
        variables = var;
        does_own_variables = false;
        functional->set_variables(var);
      }
      virtual void destroy_variables()
      {
        functional->destroy_variables();
        variables = NULL;
        does_own_variables = false;
      }
      virtual t_Container* get_variables() const
        { return  variables; }

      virtual void resize(const types::t_unsigned _nb) 
      {
        functional->resize(_nb);
        variables = functional->get_variables();
        does_own_variables = false;
      }

      virtual bool init()
       { return functional->init(); }

      virtual bool is_taboo() const
       { return functional->is_taboo(); }
  };
  template<typename T_FUNC1, typename T_FUNC2>
  class Binary_Functor : public Base<typename T_FUNC1::t_Type, typename T_FUNC2::t_Container>
  {
    public:
      typedef T_FUNC1 t_Functional1;
      typedef T_FUNC2 t_Functional2;
      typedef typename t_Functional1 :: t_Type t_Type;
      typedef typename t_Functional1 :: t_Container t_Container;
      using Base<t_Type, t_Container> :: size;

    protected:
      using Base<t_Type, t_Container> :: variables;
      using Base<t_Type, t_Container> :: does_own_variables;

    protected:
      t_Functional1 *functional1;
      t_Functional2 *functional2;

    public:     // t_Functional1 and t_Functional2 should exist, so doesn't call
      explicit  // constructors and destructors
        Binary_Functor   (const t_Functional1 *_functional1,
                          const t_Functional2 *_functional2)
                       : functional1(_functional1), functional2(_functional2)
        { 
          variables = functional1->get_variables();
          does_own_variables = false;
          functional2->set_variables(variables); 
        };
      Binary_Functor(const Binary_Functor<t_Functional1, t_Functional2> &_f) 
      {
        functional1 = _f.functional1; functional2 = _f.functional2;
        variables = functional1->get_variables();
        does_own_variables = false;
        functional2->set_variables(variables); 
      };
      Binary_Functor() : functional1(NULL), functional2(NULL)
        { variables = NULL; does_own_variables = false; } 
                   
      virtual ~Binary_Functor() {};
      
      void set_functional1( t_Functional1 * const _functional1, bool set_variables_to_obj = true )
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
      void set_functional2( t_Functional2 *const _functional2 )
      { 
        functional2 = _functional2; 
        functional2->set_variables(variables);
      }
      t_Functional1* get_functional1() const { return functional1; };
      t_Functional2* get_functional2() const { return functional2; };
    
      // required minimizer::?? behavior
    public:
      virtual void set_variables( t_Container* const var ) 
      {
        variables = var;
        does_own_variables = false;
        functional1->set_variables(var);
        functional2->set_variables(var);
      }
      virtual void destroy_variables()
      {
        if (does_own_variables and variables)
          delete variables;
        variables = NULL; does_own_variables = false;
        functional1->set_variables(NULL);
        functional2->set_variables(NULL);
      }
      virtual t_Container* get_variables() const
        { return  variables; }

      virtual void resize(const types::t_unsigned _nb) 
      {
        functional1->resize(_nb);
        variables = functional1->get_variables();
        does_own_variables = false;
        functional2->set_variables(variables);
      }

      virtual bool init()
       { return functional1->init() and functional2->init(); }

      virtual bool is_taboo() const
       { return functional1->is_taboo() or functional2->is_taboo(); }
  };

  template<typename t_Functional1, typename t_Functional2>
  class Minus : public Binary_Functor<t_Functional1, t_Functional2>
  {
    public:
      typedef typename Binary_Functor<t_Functional1, t_Functional2> :: t_Type t_Type;
      typedef typename Binary_Functor<t_Functional1, t_Functional2> :: t_Container t_Container;
      using Binary_Functor<t_Functional1, t_Functional2> :: size;

    protected:
      using Binary_Functor<t_Functional1, t_Functional2> :: functional1;
      using Binary_Functor<t_Functional1, t_Functional2> :: functional2;
      using Binary_Functor<t_Functional1, t_Functional2> :: variables;
      using Binary_Functor<t_Functional1, t_Functional2> :: does_own_variables;


    public:     // t_Functional1 and t_Functional2 should exist, so doesn't call
      explicit  // constructors and destructors
        Minus(const t_Functional1 *_functional1, const t_Functional2 *_functional2)
        { 
          functional1 = _functional1; functional2 = _functional2;
          variables = functional1->get_variables(); does_own_variables = false;
          functional2->set_variables(variables); 
        };
      Minus() : Binary_Functor<t_Functional1, t_Functional2>() {};
      virtual ~Minus() {};
      
      // required minimizer::?? behavior
    public:
      virtual t_Type evaluate() 
        { return functional1->evaluate() - functional2->evaluate(); }
      virtual void evaluate_gradient( t_Type* const _i_grad) 
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
      virtual t_Type evaluate_one_gradient( types::t_unsigned _pos) 
      {
        return   functional1->evaluate_one_gradient(_pos)
               + functional2->evaluate_one_gradient(_pos);
      }
      virtual t_Type evaluate_with_gradient( t_Type* const _i_grad)
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
  };

  template<typename t_Functional1, typename t_Functional2>
  class Plus : public Binary_Functor<t_Functional1, t_Functional2>
  {
    public:
      typedef typename Binary_Functor<t_Functional1, t_Functional2> :: t_Type t_Type;
      typedef typename Binary_Functor<t_Functional1, t_Functional2> :: t_Container t_Container;
      using Binary_Functor<t_Functional1, t_Functional2> :: size;

    protected:
      using Binary_Functor<t_Functional1, t_Functional2> :: functional1;
      using Binary_Functor<t_Functional1, t_Functional2> :: functional2;
      using Binary_Functor<t_Functional1, t_Functional2> :: variables;
      using Binary_Functor<t_Functional1, t_Functional2> :: does_own_variables;

    public:     // t_Functional1 and t_Functional2 should exist, so doesn't call
      explicit  // constructors and destructors
        Plus(const t_Functional1 *_functional1, const t_Functional2 *_functional2)
        { 
          functional1 = _functional1;  functional2 = _functional2; 
          variables = functional1->get_variables(); does_own_variables = false;
          functional2->set_variables(variables); 
        };
      Plus() : Binary_Functor<t_Functional1, t_Functional2>() {};
      virtual ~Plus(){};
      
      // required minimizer:: behavior
    public:
      virtual t_Type evaluate() 
        { return (functional1->evaluate() + functional2->evaluate()); }
      virtual void evaluate_gradient( t_Type* const _i_grad) 
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
      virtual t_Type evaluate_one_gradient( types::t_unsigned _pos) 
      {
        return   functional1->evaluate_one_gradient(_pos)
               - functional2->evaluate_one_gradient(_pos);
      }
      virtual t_Type evaluate_with_gradient( t_Type* const _i_grad ) 
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
  };

  template<typename T_FUNCTIONAL>
  class Quadratic : public Unary_Functor< T_FUNCTIONAL >
  {
    public:
      typedef T_FUNCTIONAL t_Functional;
      typedef typename t_Functional :: t_Type t_Type;
      typedef typename t_Functional :: t_Container t_Container;

    protected:
      using Unary_Functor<t_Functional> :: functional;
      using Base<t_Type, t_Container> :: size;
      t_Type center;

    public:     // t_Functional1 and t_Functional2 should exist, so doesn't call
      explicit  // constructors and destructors
        Quadratic   (const t_Functional *_func, const t_Type &_center = t_Type(0) )
                  : Unary_Functor<t_Functional>(_func),
                    center(_center) {}
      Quadratic   ( const t_Type &_center = t_Type(0) )
                : function::Unary_Functor<t_Functional>(), center(_center) {};
      virtual ~Quadratic(){};
      
      // required minimizer::? behavior
    public:
      virtual t_Type evaluate() 
        { t_Type value = functional->evaluate() - center; return value * value; }
      virtual t_Type evaluate_with_gradient( t_Type* const _i_grad) 
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
      virtual void evaluate_gradient( t_Type* const _i_grad ) 
        { evaluate_with_gradient(_i_grad); }
      virtual t_Type evaluate_one_gradient( types::t_unsigned _pos) 
      {
        t_Type value = functional->evaluate() - center;
        return  2.0 * functional->evaluate_one_gradient(_pos) * value;
      }
  };
} // namespace function
#endif // _OPT_FUNCTORS_H

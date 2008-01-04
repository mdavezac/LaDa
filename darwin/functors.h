//
//  Version: $Id$
//
#ifndef _DARWIN_FUNCTORS_H_
#define _DARWIN_FUNCTORS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <eo/eoOp.h>
#include <eo/eoContinue.h>
#include <eo/utils/eoRNG.h>
#include <eo/utils/eoState.h>

#include <opt/types.h>
#include <lamarck/structure.h>

#include "gatraits.h"

namespace GA 
{
  //! \brief Defines a constant unary EO functor
  //! \details Can be called by const member %functions. Can be used with
  //!          eostates type storage.
  template <class A1, class R>
 class const_eoUF : public eoFunctorBase, public std::unary_function<A1, R>
 {
   public:
     //! Type of the argument
     typedef A1 t_Argument;
     //! Type of the returned value
     typedef R  t_Return;

   public:
     //! The virtual destructor
     virtual ~const_eoUF() {}
   
     //! The functor member itself
     virtual t_Return operator()( t_Argument ) const = 0;
   
     //! The functor category
     static eoFunctorBase::unary_function_tag functor_category()
       { return eoFunctorBase::unary_function_tag(); }
 };

  //! \brief generic class for converting GA::Evaluator member %function to
  //!        binary operators
  template<class T_GATRAITS >
  class mem_monop_t : public eoMonOp<typename T_GATRAITS::t_Individual > 
  {
    public:
      //! All relevant %GA traits.
      typedef T_GATRAITS t_GATraits;
    protected:
      //! The evaluator type
      typedef typename t_GATraits :: t_Evaluator t_Evaluator; 
      //! The type of the individual
      typedef typename t_GATraits :: t_Individual t_Individual; 
    public:
      //! Pointer to a member %function
      typedef bool ( t_Evaluator :: *t_Function )( t_Individual &);

    private:
      //! Reference to a GA::Evaluator derived instance
      t_Evaluator &class_obj;
      //! The member %function for which to act as a functor
      t_Function class_func;
      //! A string characterizing the new functor
      std::string class_name;

    public:
      //! Constructor and Initializer
      explicit
        mem_monop_t   ( t_Evaluator &_co, t_Function _func, const std::string &_cn )
                    : class_obj(_co), class_func(_func), class_name(_cn) {};

      //! Sets the string characterizing this functor
      void set_className( std::string &_cn) { class_name = _cn; }
      //! Returns the string characterizing this functor
      virtual std::string className() const { return class_name; }

      //! The functor itself
      bool operator() (t_Individual &_object) 
        {  return ( (class_obj.*class_func) )( _object); }

  }; 
  //! generic class for converting member %function to zero operators
  template<class T_CLASS>
  class mem_zerop_t : public eoF<bool>
  {
    public:
      //! The type of the class for which to construct a functor
      typedef T_CLASS t_Class; 
      //! Pointer to a member %function
      typedef bool ( t_Class::*t_Function )();

    private:
      //! Reference to an instance of the calss
      t_Class &class_obj;
      //! The member %function for which to act as a functor
      t_Function class_func;
      //! A string characterizing the new functor
      std::string class_name;

    public:
      //! Constructor and Initializer
      explicit
        mem_zerop_t   ( t_Class &_co, t_Function _func, const std::string &_cn )
                    : class_obj(_co), class_func(_func), class_name(_cn) {};

      //! Sets the string characterizing this functor
      void set_className( std::string &_cn) { class_name = _cn; }
      //! Returns the string characterizing this functor
      virtual std::string className() const { return class_name; }

      //! The functor itself
      bool operator() () 
        {  return ( (class_obj.*class_func) )(); }

  }; 
  
  //! \brief generic class for converting %GA::Evaluator member %function with
  //!        an argument to binary genetic operators
  template<class T_EVALUATOR, class T_ARG>
  class mem_bingenop_arg_t : public eoGenOp<typename T_EVALUATOR :: t_Individual >
  {
    public:
      //! The evaluator type
      typedef T_EVALUATOR t_Evaluator;
      //! The type of the argument
      typedef T_ARG t_Arg; 
    protected:
      //! The type of the individual
      typedef typename t_Evaluator :: t_Individual t_Individual; 
    public:
      //! Pointer to a member %function
      typedef bool ( t_Evaluator::*t_Function )(t_Individual &, const t_Individual&, t_Arg);

    private:
      //! Reference to a GA::Evaluator derived instance
      t_Evaluator &class_obj;
      //! The argument with which to call the member %function
      t_Arg arg;
      //! The member %function for which to act as a functor
      t_Function class_func;
      //! A string characterizing the new functor
      std::string class_name;

    public:
      //! Constructor and Initializer
      explicit
        mem_bingenop_arg_t   ( t_Evaluator &_co, t_Function _func,
                        const std::string &_cn, t_Arg _arg )
                    : class_obj(_co), arg(_arg), class_func(_func), class_name(_cn) {};

      //! Sets the string characterizing this functor
      void set_className( std::string &_cn) { class_name = _cn; }
      //! Returns the string characterizing this functor
      virtual std::string className() const { return class_name; }

      //! The number of created individuals
      unsigned max_production(void) { return 1; } 

      //! Makes this class a genetic operator
      void apply(eoPopulator<t_Individual>& _pop)
      {
        t_Individual& a = *_pop;
        const t_Individual& b = _pop.select();
  
        if ( (class_obj.*class_func)(a, b, arg) )
          a.invalidate();
      }
  }; 

  //! \brief generic class for converting GA::Evaluator member %function to binary
  //!        genetic operators
  template<class T_EVALUATOR>
  class mem_bingenop_t : public eoGenOp<typename T_EVALUATOR :: t_Individual> 
  {
    public:
      //! The evaluator type
      typedef T_EVALUATOR t_Evaluator;
    protected:
      //! The type of the individual
      typedef typename t_Evaluator::t_Individual t_Individual; 
    public:
      //! Pointer to a member %function
      typedef bool ( t_Evaluator::*t_Function )(t_Individual &, const t_Individual&);

    private:
      //! Reference to a GA::Evaluator derived instance
      t_Evaluator &class_obj;
      //! The member %function for which to act as a functor
      t_Function class_func;
      //! A string characterizing the new functor
      std::string class_name;

    public:
      //! Reference to a GA::Evaluator derived instance
      explicit
        mem_bingenop_t   ( t_Evaluator &_co, t_Function _func,
                           const std::string &_cn )
                    : class_obj(_co), class_func(_func), class_name(_cn) {};

      //! Sets the string characterizing this functor
      void set_className( std::string &_cn) { class_name = _cn; }
      //! Returns the string characterizing this functor
      virtual std::string className() const { return class_name; }

      //! The number of created individuals
      unsigned max_production(void) { return 1; } 

      //! Makes this class a genetic operator
      void apply(eoPopulator<t_Individual>& _pop)
      {
        t_Individual& a = *_pop;
        const t_Individual& b = _pop.select();
  
        if ( (class_obj.*class_func)(a, b) )
          a.invalidate();
      }
  }; 

  //! \brief Helper %function for wrapping an GA::Evaluator member %function into a
  //!        genetic operator.
  template< class T_EVALUATOR, class T_ARGS >
    mem_bingenop_arg_t<T_EVALUATOR, T_ARGS>*
        new_genop( T_EVALUATOR &_eval, 
                   typename mem_bingenop_arg_t<T_EVALUATOR, T_ARGS> :: t_Function _func,
                   const std::string  &_str, T_ARGS _arg )
    {
      return new mem_bingenop_arg_t<T_EVALUATOR, T_ARGS>( _eval, _func, _str, _arg );
    }
  //! \brief Helper %function for wrapping an GA::Evaluator member %function into a
  //!        genetic operator.
  template< class T_EVALUATOR >
    mem_bingenop_t<T_EVALUATOR>* 
      new_genop ( T_EVALUATOR &_eval,
                  typename mem_bingenop_t<T_EVALUATOR> :: t_Function _func,
                  const std::string  &_str )
    {
      return new mem_bingenop_t<T_EVALUATOR>( _eval, _func, _str );
    }


  //! a dummy unary operator which does nothing 
  template< class T_OBJECT >
  class DummyOp : public eoMonOp<T_OBJECT>
  {
    protected:
      //! Type of the argument
      typedef T_OBJECT t_Object;

    public:
      //! Constructor
      DummyOp() {};
      //! Destructor
      ~DummyOp() {}

      //! Does nothing!
      bool operator() ( t_Object &_obj1 )
      { return false; } // do nothing!
  };

  //! \cond
  // Historycal shit ?
  template< class T_GATRAITS >
  class Continuator : public eoContinue< typename T_GATRAITS::t_Individual >
  {
    public:
      //! All relevant %GA traits
      typedef T_GATRAITS t_GATraits;
    protected:
      //! Type of the individual
      typedef typename t_GATraits::t_Individual t_Individual;
      //! Type of the population
      typedef typename t_GATraits::t_Population t_Population;

    protected:
      //! Reference
      eoF<bool> &op;

    public:
      Continuator( eoF<bool> &_op ) : op(_op) {};
      ~Continuator() {}

      bool operator()(const t_Population &_pop )
      {
        return op();
      }
  };
  //! \endcond


}
#endif // _DARWIN_FUNCTORS_H_

//
//  Version: $Id$
//
#ifndef _COMPOSE_FUNCTORS_H
#define _COMPOSE_FUNCTORS_H
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <functional>
namespace opt 
{
  //! \brief Composes three functors \f$M(x,y)\f$, \f$g_0(x)\f$, and
  //!        \f$g_1(y)\f$ to create the functor \f$M(\ g_0(x),\ g_1(y)\ )\f$.
  //! \details In practice, you should use the function compose_binary() to
  //!          obtain an object of this kind. The functors explicitely takes
  //           references as arguments.
  template <class _Operation1, class _Operation2, class _Operation3>
    class binary_compose_binary : public std::binary_function<typename _Operation2::argument_type,
			                                      typename _Operation3::argument_type,
			                                      typename _Operation1::result_type>
    {
    protected:
      //! Binary Functor 
      _Operation1 _M_fn1;
      //! Unary Functor, first argument to binary.
      _Operation2 _M_fn2;
      //! Unary Functor, second argument to binary.
      _Operation3 _M_fn3;
      
    public:
      //! Constructor
      binary_compose_binary(const _Operation1& __x, const _Operation2& __y,
		            const _Operation3& __z)
      : _M_fn1(__x), _M_fn2(__y), _M_fn3(__z) { }

      //! Functor itself
      typename _Operation1::result_type
      operator()(const typename _Operation2::argument_type& __x,
                 const typename _Operation3::argument_type& __y) const
      { return _M_fn1(_M_fn2(__x), _M_fn3(__y)); }
    };


  //! \brief Composes three functors \f$M(x,y)\f$, \f$g_0(x)\f$, and
  //!        \f$g_1(y)\f$ to create the functor \f$M(\ g_0(x),\ g_1(y)\ )\f$.
  //! \details In practice, you should use the function compose_binary() to
  //!          obtain an object of this kind.
  template <class _Operation1, class _Operation2, class _Operation3>
  class ref_binary_compose_binary
    : public std::binary_function<typename _Operation3::argument_type,
                                  typename _Operation2::argument_type,
                                  typename _Operation1::result_type>
  {
    protected:
      //! Binary Functor 
      _Operation1 _M_fn1;
      //! Unary Functor, first argument to binary.
      _Operation2 _M_fn2;
      //! Unary Functor, second argument to binary.
      _Operation3 _M_fn3;
      
    public:
      //! Constructor
      ref_binary_compose_binary   (const _Operation1& __x, const _Operation2& __y,
                                   const _Operation3& __z)
                                : _M_fn1(__x), _M_fn2(__y), _M_fn3(__z) { }

      //! Functor itself
      typename _Operation1::result_type
      operator()(const typename _Operation2::argument_type __x,
                 const typename _Operation3::argument_type __y) const
      { return _M_fn1(_M_fn2(__x), _M_fn3(__y)); }
  };

  //! A helper functions for composing functor \f$M(\ g_0(x),\ g_1(y)\ )\f$.
  //! \param __fn1 binary functor \f$M(x,y)\f$
  //! \param __fn2 unary functor \f$g_0(x)\f$
  //! \param __fn3 unary functor \f$g_1(y)\f$
  template <class _Operation1, class _Operation2, class _Operation3>
    inline ref_binary_compose_binary<_Operation1, _Operation2, _Operation3>
      ref_compose2(const _Operation1& __fn1, const _Operation2& __fn2,
                   const _Operation3& __fn3)
      { 
        return ref_binary_compose_binary<_Operation1, _Operation2, _Operation3>(__fn1, __fn2, __fn3); 
      }
  
  //! \brief Creates a composed functor \f$ f( g( x ) ) \f$, where the functor
  //!        takes explicitely a reference,
  template <class _Operation1, class _Operation2>
  class ref_unary_compose
    : public std::unary_function<typename _Operation1::argument_type,
                                 typename _Operation2::result_type>
  {
    protected:
      //! Outer unary functor.
      _Operation1 _M_fn1;
      //! Inner unary functor.
      _Operation2 _M_fn2;
      
    public:
      //! \brief Constructor
      //! \param __x Outer functor.
      //! \param __y Inner functor.
      ref_unary_compose   (const _Operation1& __x, const _Operation2& __y )
                        : _M_fn1(__x), _M_fn2(__y) {}

      //! Functor itself. Explicitely takes a reference.
      typename _Operation1::result_type
      operator()(const typename _Operation2::argument_type __x) const
      { return _M_fn1(_M_fn2(__x)); }
  };

  //! \brief Helper function  for composing unary functors which take
  //!        references.
  //! \param __fn1 Outer functor.
  //! \param __fn2 Inner functor.
  template <class _Operation1, class _Operation2>
    inline ref_unary_compose<_Operation1, _Operation2>
      ref_compose1(const _Operation1& __fn1, const _Operation2& __fn2 )
      { 
        return ref_unary_compose<_Operation1, _Operation2>(__fn1, __fn2); 
      }

  //! \brief Creates a unary functor from a binary %function by binding the
  //!        second argument.
  //! \details The binary %function explicitely takes references for both
  //!          arguments. As a result, so fdoes this functor.
  template <class _Arg1, class _Arg2, class _Ret>
  class const_binder2nd_ref  : public std::unary_function<_Arg1, _Ret>
  {
    //! Type of the pointer to a binary function.
    typedef _Ret (*_Operation)(_Arg1&, _Arg2& );
    protected:
      //! Pointer to the binary function.
      _Operation _M_fn;
      //! Bound reference and second argument to the binary %function. 
      _Arg2 &_boundArg;
      
    public:
      //! \brief Constructor.
      //! \param __x Adress of the binary %function.
      //! \param __y Reference to the second (bound) argument.
      const_binder2nd_ref   (const _Operation  __x, _Arg2& __y )
                          : _M_fn(__x), _boundArg(__y) {}

      //! The functor itself.
      _Ret operator()(const _Arg1& __x) const
        { return _M_fn( __x, _boundArg ); }
  };
  //! \brief Helper function for binding to a reference the second argument of a
  //!        binary %function.
  //! \param __fn Adress of the binary %function.
  //! \param __x Reference to the second (bound) argument.
  template <class _Arg1, class _Arg2, class _Ret>
    inline const_binder2nd_ref<_Arg1, _Arg2, _Ret>
    bind2nd(_Ret (*__fn)(_Arg1&, _Arg2&), _Arg2 &__x)
      { return const_binder2nd_ref<_Arg1, _Arg2, _Ret>(__fn, __x); }



  //! \brief Creates a binary functor for ternary %function by binding the
  //!        third argument to a reference.
  template <class _Arg1, class _Arg2, class _Arg3, class _Ret>
  class binder3rd_ref  : public std::binary_function<_Arg1, _Arg2, _Ret>
  {
    //! Type of the pointer to a ternary function.
    typedef _Ret (*_Operation)(_Arg1&, _Arg2&, _Arg3 );
    protected:
      //! Pointer to the ternary function.
      _Operation _M_fn;
      //! Bound reference and third argument to the ternary %function. 
      _Arg3 &_boundArg;
      
    public:
      //! \brief Constructor
      //! \param __x Adress of the ternary %function.
      //! \param __y Reference to the third (bound) argument.
      binder3rd_ref   (const _Operation  __x, const _Arg3& __y )
                    : _M_fn(__x), _boundArg(__y) {}

      //! The functor itself. Takes two references, explicitely.
      _Ret operator()(const _Arg1& __x, const _Arg2& __y) const
        { return _M_fn( __x, __y, _boundArg ); }
  };
  //! \brief Helper function for binding to a reference the third argument of a
  //!        ternary %function.
  //! \param __fn Adress of the ternary %function.
  //! \param __x Reference to the third (bound) argument.
  template <class _Arg1, class _Arg2, class _Arg3, class _Ret, class _Tp>
    inline binder3rd_ref<_Arg1, _Arg2, _Arg3, _Ret>
      bind3rd(_Ret (*__fn)(_Arg1&, _Arg2&, _Arg3& ), const _Tp __x)
        { return binder3rd_ref<_Arg1, _Arg2, _Arg3, _Ret>(__fn, _Arg3(__x)); }
}

#endif

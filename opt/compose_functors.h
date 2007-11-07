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
  //!          obtain an object of this kind.
  template <class _Operation1, class _Operation2, class _Operation3>
    class binary_compose_binary
    : public binary_function<typename _Operation2::argument_type,
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

  //! A helper functions for composing functor \f$M(\ g_0(x),\ g_1(y)\ )\f$.
  //! \param __fn1 binary functor \f$M(x,y)\f$
  //! \param __fn2 unary functor \f$g_0(x)\f$
  //! \param __fn3 unary functor \f$g_1(y)\f$
  template <class _Operation1, class _Operation2, class _Operation3>
    inline binary_compose_binary<_Operation1, _Operation2, _Operation3>
    compose_binary(const _Operation1& __fn1, const _Operation2& __fn2,
	           const _Operation3& __fn3)
    { return binary_compose_binary<_Operation1, _Operation2, _Operation3>
	(__fn1, __fn2, __fn3); }



  template <class _Operation1, class _Operation2, class _Operation3>
  class ref_binary_compose_binary
    : public std::binary_function<typename _Operation3::argument_type,
                                  typename _Operation2::argument_type,
                                  typename _Operation1::result_type>
  {
    protected:
      _Operation1 _M_fn1;
      _Operation2 _M_fn2;
      _Operation3 _M_fn3;
      
    public:
      ref_binary_compose_binary   (const _Operation1& __x, const _Operation2& __y,
                                   const _Operation3& __z)
                                : _M_fn1(__x), _M_fn2(__y), _M_fn3(__z) { }

      typename _Operation1::result_type
      operator()(const typename _Operation2::argument_type __x,
                 const typename _Operation3::argument_type __y) const
      { return _M_fn1(_M_fn2(__x), _M_fn3(__y)); }
  };

  template <class _Operation1, class _Operation2, class _Operation3>
    inline ref_binary_compose_binary<_Operation1, _Operation2, _Operation3>
      ref_compose2(const _Operation1& __fn1, const _Operation2& __fn2,
                   const _Operation3& __fn3)
      { 
        return ref_binary_compose_binary<_Operation1, _Operation2, _Operation3>(__fn1, __fn2, __fn3); 
      }
  
  template <class _Operation1, class _Operation2>
  class ref_unary_compose
    : public std::unary_function<typename _Operation1::argument_type,
                                 typename _Operation2::result_type>
  {
    protected:
      _Operation1 _M_fn1;
      _Operation2 _M_fn2;
      
    public:
      ref_unary_compose   (const _Operation1& __x, const _Operation2& __y )
                        : _M_fn1(__x), _M_fn2(__y) {}

      typename _Operation1::result_type
      operator()(const typename _Operation2::argument_type __x) const
      { return _M_fn1(_M_fn2(__x)); }
  };

  template <class _Operation1, class _Operation2>
    inline ref_unary_compose<_Operation1, _Operation2>
      ref_compose1(const _Operation1& __fn1, const _Operation2& __fn2 )
      { 
        return ref_unary_compose<_Operation1, _Operation2>(__fn1, __fn2); 
      }

  template <class _Arg1, class _Arg2, class _Ret>
  class const_binder2nd_ref  : public std::unary_function<_Arg1, _Ret>
  {
    typedef _Ret (*_Operation)(_Arg1&, _Arg2& );
    protected:
      _Operation _M_fn;
      _Arg2 &_boundArg;
      
    public:
      const_binder2nd_ref   (const _Operation  __x, _Arg2& __y )
                          : _M_fn(__x), _boundArg(__y) {}

      _Ret operator()(const _Arg1& __x) const
        { return _M_fn( __x, _boundArg ); }
  };
  template <class _Arg1, class _Arg2, class _Ret>
    inline const_binder2nd_ref<_Arg1, _Arg2, _Ret>
    bind2nd(_Ret (*__fn)(_Arg1&, _Arg2&), _Arg2 &__x)
    {
      return const_binder2nd_ref<_Arg1, _Arg2, _Ret>(__fn, __x);
    }
  template <class _Arg1, class _Arg2, class _Arg3, class _Ret>
  class binder3rd_ref  : public std::binary_function<_Arg1, _Arg2, _Ret>
  {
    typedef _Ret (*_Operation)(_Arg1&, _Arg2&, _Arg3 );
    protected:
      _Operation _M_fn;
      _Arg3 _boundArg;
      
    public:
      binder3rd_ref   (const _Operation  __x, const _Arg3& __y )
                    : _M_fn(__x), _boundArg(__y) {}

      _Ret operator()(const _Arg1& __x, const _Arg2& __y) const
        { return _M_fn( __x, __y, _boundArg ); }
  };
  template <class _Arg1, class _Arg2, class _Arg3, class _Ret, class _Tp>
    inline binder3rd_ref<_Arg1, _Arg2, _Arg3, _Ret>
    bind3rd(_Ret (*__fn)(_Arg1&, _Arg2&, _Arg3 ), const _Tp __x)
    {
      return binder3rd_ref<_Arg1, _Arg2, _Arg3, _Ret>(__fn, _Arg3(__x));
    }
}

#endif

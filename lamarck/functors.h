//
//  Version: $Id$
//
#ifndef _FUNCTORS_H_
#define _FUNCTORS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<functional>

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

#endif

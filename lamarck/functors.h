#ifndef _FUNCTORS_H_
#define _FUNCTORS_H_

#include<functional>

  template <class _Operation1, class _Operation2, class _Operation3>
    class binary_compose_binary
    : public binary_function<typename _Operation2::argument_type,
			     typename _Operation3::argument_type,
			     typename _Operation1::result_type>
    {
    protected:
      _Operation1 _M_fn1;
      _Operation2 _M_fn2;
      _Operation3 _M_fn3;
      
    public:
      binary_compose_binary(const _Operation1& __x, const _Operation2& __y,
		            const _Operation3& __z)
      : _M_fn1(__x), _M_fn2(__y), _M_fn3(__z) { }

      typename _Operation1::result_type
      operator()(const typename _Operation2::argument_type& __x,
                 const typename _Operation3::argument_type& __y) const
      { return _M_fn1(_M_fn2(__x), _M_fn3(__y)); }
    };

  /// An \link SGIextensions SGI extension \endlink.
  template <class _Operation1, class _Operation2, class _Operation3>
    inline binary_compose_binary<_Operation1, _Operation2, _Operation3>
    compose_binary(const _Operation1& __fn1, const _Operation2& __fn2,
	           const _Operation3& __fn3)
    { return binary_compose_binary<_Operation1, _Operation2, _Operation3>
	(__fn1, __fn2, __fn3); }

#endif

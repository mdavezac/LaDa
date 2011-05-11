//
//  Version: $Id: protect.h 1189 2009-06-17 02:19:52Z davezac $
//

#ifndef _LADA_LOADNSAVE_TRANSFORM_WRAP_H_
#define _LADA_LOADNSAVE_TRANSFORM_WRAP_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/type_traits/add_const.hpp>
#include <boost/type_traits/add_reference.hpp>
#include <boost/type_traits/remove_reference.hpp>
#include <boost/proto/transform/arg.hpp>
#include <boost/fusion/include/at_c.hpp>

namespace LaDa 
{
  namespace load_n_save
  {
    namespace transform
    {
      namespace proto = boost::proto;
      namespace fusion = boost::fusion;
      
      //! A primitive protecting an object with a const view.
      struct protect : proto::callable
      {
        //! Resulting type.
        template< class T_SIGNA > struct result;
        //! Specialization for calls.
        template< class T_THIS, class T_EXPR >
          struct result< T_THIS(const T_EXPR&) >
          {
            //! Returns a view.
            typedef typename fusion::single_view<const T_EXPR> const type;
          };
        //! Specialization for calls.
        template< class T_THIS, class T_EXPR >
          struct result< T_THIS(T_EXPR) >
          {
            //! Returns a view.
            typedef typename fusion::single_view<const T_EXPR> const type;
          };
        //! Specialization for calls.
        template< class T_THIS, class T_EXPR, int N >
          struct result< T_THIS(T_EXPR const(&)[N]) >
          {
            //! Returns a view.
            typedef typename fusion::single_view<T_EXPR const(&)[N]> const type;
          };
        //! Specialization for calls.
        template< class T_THIS, class T_EXPR, int N >
          struct result< T_THIS(T_EXPR (&)[N]) >
          {
            //! Returns a view.
            typedef typename fusion::single_view<T_EXPR const(&)[N]> const type;
          };
          
        //! implementation of the transform.
        template<typename T_EXPR >
          typename result<protect(const T_EXPR&)> :: type 
            operator ()( const T_EXPR& _expr ) const
            {
              return typename result<protect(const T_EXPR&)> :: type( _expr );
            }
        //! implementation of the transform.
        template<typename T_EXPR, int N >
          typename result<protect(T_EXPR const(&)[N])> :: type 
            operator ()( T_EXPR const(&_expr)[N] ) const
            {
              return typename result<protect(const T_EXPR(&)[N])> :: type( _expr );
            }
      };

      //! A primitive protecting a reference with a const view.
      struct protect_ref : proto::callable
      {
        //! Resulting type.
        template< class T_SIGNA > struct result;
        //! Specialization for calls.
        template< class T_THIS, class T_EXPR >
          struct result< T_THIS(T_EXPR) >
          {
            //! Returns a view.
            typedef typename fusion::single_view<T_EXPR&> const type;
          };
          
        //! implementation of the transform.
        template<typename T_EXPR >
          typename result<protect_ref(T_EXPR&)> :: type 
            operator ()( T_EXPR& _expr ) const
            {
              return boost::fusion::single_view< T_EXPR& >( _expr );
            }
      };

      //! A primitive protecting a const reference with a const view.
      struct protect_cref : proto::callable
      {
        //! Resulting type.
        template< class T_SIGNA > struct result;
        //! Specialization for calls.
        template< class T_THIS, class T_EXPR >
          struct result< T_THIS(T_EXPR) >
          {
            //! Returns a view.
            typedef typename fusion::single_view<T_EXPR const &> const type;
          };
          
        //! implementation of the transform.
        template<typename T_EXPR >
          typename result<protect_ref(T_EXPR)> :: type 
            operator ()( const T_EXPR& _expr ) const
            {
              return boost::fusion::single_view< T_EXPR const & >( _expr );
            }
      };

    } // namespace transform
  } // namespace load_n_save

} // namespace LaDa

#endif

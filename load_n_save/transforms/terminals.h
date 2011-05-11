//
//  Version: $Id: terminals.h 1212 2009-07-04 04:36:19Z davezac $
//

#ifndef _LADA_LOADNSAVE_TRANSFORM_TERMINALS_H_
#define _LADA_LOADNSAVE_TRANSFORM_TERMINALS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/proto/transform/arg.hpp>
#include <boost/fusion/include/at_c.hpp>

namespace LaDa 
{
  namespace load_n_save
  {
    namespace transform
    {
      namespace proto = boost::proto;
      
      //! A primitive transform returning the name (member variable) of a terminal.
      struct name : proto::callable
      {
        //! Resulting type is a const char*.
        typedef const char * result_type;
        //! implementation of the transform.
        template<typename T_EXPR >
          result_type operator ()( const T_EXPR& _expr ) const
            { return result_type( proto::child( _expr ).value ); }
      };

      //! A primitive transform returning the name (member variable) of a terminal.
      struct tag : proto::callable
      {
        //! Resulting type is a const char*.
        typedef size_t result_type;
        //! implementation of the transform.
        template<typename T_EXPR >
          result_type operator ()( const T_EXPR& _expr ) const
            { return result_type( proto::child( _expr ).tag ); }
      };


      //! A primitive transform returning the callable const ref.
      struct var : proto::callable
      {
        //! Resulting type.
        template< class T > struct result;
        //! Specialization for calls.
        template< class T_THIS, class T_EXPR >
          struct result<T_THIS(T_EXPR)> 
            : public boost::fusion::result_of::at_c<T_EXPR, 0> {};
        //! Specialization for calls.
        template< class T_THIS, class T_EXPR >
          struct result<T_THIS(T_EXPR&)> 
            : public boost::fusion::result_of::at_c<T_EXPR, 0> {};
        //! Specialization for vars.
        template< class T_THIS, class T_EXPR >
          struct result<T_THIS(const T_EXPR&)> 
            : public boost::fusion::result_of::at_c<T_EXPR, 0> {};
          
        //! implementation of the transform.
        template<typename T_EXPR >
          typename result<var(const T_EXPR&)> :: type 
            operator ()( T_EXPR& _expr ) const
              { return boost::fusion::at_c<0>(_expr); }
      };
      //! A primitive transform returning the callable const ref.
      struct call : proto::callable
      {
        //! Resulting type.
        template< class T > struct result;
        //! Specialization for calls.
        template< class T_THIS, class T_EXPR >
          struct result<T_THIS(T_EXPR&)> 
            : public boost::add_const
              < 
                typename boost::fusion::result_of::at_c<const T_EXPR, 0> :: type 
              > {};
        //! Specialization for vars.
        template< class T_THIS, class T_EXPR >
          struct result<T_THIS(const T_EXPR&)> 
            : public boost::add_const
              < 
                typename boost::fusion::result_of::at_c<const T_EXPR, 0> :: type 
              > {};
          
        //! implementation of the transform.
        template<typename T_EXPR >
          typename result<call(const T_EXPR&)> :: type 
            operator ()( const T_EXPR& _expr ) const
              { return boost::fusion::at_c<0>(_expr); }
      };
    } // namespace transform
  } // namespace load_n_save

} // namespace LaDa

#endif

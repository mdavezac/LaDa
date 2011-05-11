//
//  Version: $Id: join.h 1128 2009-05-20 04:38:42Z Mayeul $
//

#ifndef _LADA_LOADNSAVE_TRANSFORM_JOIN_H_
#define _LADA_LOADNSAVE_TRANSFORM_JOIN_H_

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
      
      //! A primitive transform returning the callable const ref.
      struct join : proto::callable
      {
        //! Resulting type.
        template< class T_SIGNA > struct result;
        //! Specialization for binary calls.
        template< class T_THIS, class T1, class T2 >
          struct result< T_THIS(T1, T2) >
          {
            //! Returns a view.
            typedef typename fusion::joint_view<const T1, const T2> const type;
          };
        //! Specialization for ternary calls.
        template< class T_THIS, class T1, class T2, class T3 >
          struct result< T_THIS(T1, T2, T3) >
          {
            //! Returns a view.
            typedef typename fusion::joint_view
            <
              typename result<T_THIS(const T1, const T2)> :: type,
              const T3
            > const type;
          };
        //! Specialization for ternary calls.
        template< class T_THIS, class T1, class T2, class T3, class T4 >
          struct result< T_THIS(T1, T2, T3, T4) >
          {
            //! Returns a view.
            typedef typename fusion::joint_view
            <
              typename result<T_THIS(const T1, const T2)> :: type,
              typename result<T_THIS(const T3, const T4)> :: type 
            > const type;
          };
          
        //! implementation of the transform.
        template<class T1, class T2 >
          typename result<join(T1, T2)> :: type 
            operator ()( const T1& _t1, const T2& _t2 ) const
            {
              typedef typename result<join(T1, T2)> :: type t_Result;
              return t_Result( _t1, _t2 );
            }

        //! implementation of the transform.
        template<class T1, class T2, class T3 >
          typename result<join(T1, T2, T3)> :: type 
            operator ()( const T1& _t1, const T2& _t2, const T3& _t3 ) const
            {
              typedef typename result<join(T1, T2, T3)> :: type t_Result;
              typedef typename result<join(T1, T2)> :: type t_first;
              return t_Result( t_first(_t1, _t2), _t3 );
            }

        //! implementation of the transform.
        template<class T1, class T2, class T3, class T4 >
          typename result<join(T1, T2, T3, T4)> :: type 
            operator ()( const T1& _t1, const T2& _t2,
                         const T3& _t3, const T4 &_t4 ) const
            {
              typedef typename result<join(T1, T2)> :: type t_first;
              typedef typename result<join(T3, T4)> :: type t_second;
              typedef typename result<join(t_first, t_second)> :: type t_Result;
              return t_Result( t_first(_t1, _t2), t_second(_t3, _t4) );
            }
      };
    } // namespace transform
  } // namespace load_n_save

} // namespace LaDa

#endif

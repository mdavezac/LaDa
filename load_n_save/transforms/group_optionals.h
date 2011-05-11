//
//  Version: $Id: group_optionals.h 1226 2009-07-13 06:28:01Z davezac $
//

#ifndef _LADA_LOADNSAVE_TRANSFORMS_GROUP_OPTIONALS_H_
#define _LADA_LOADNSAVE_TRANSFORMS_GROUP_OPTIONALS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/fusion/include/begin.hpp>
#include <boost/fusion/include/is_view.hpp>
#include <boost/fusion/include/single_view.hpp>
#include <boost/fusion/include/joint_view.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/fusion/include/as_vector.hpp>
#include <boost/proto/transform/arg.hpp>
#include <boost/type_traits/remove_all_extents.hpp>

#include "join.h"
#include "protect.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace transform
    {

      template< class T >
        struct OptionalGroups
        {
          typedef typename boost::remove_all_extents<T>::type const type;
          OptionalGroups( type _t )  : view(_t) {};
          type view;
        };

      struct grouped_options : boost::proto::callable
      {
        private:
          template< class T >
            struct get_view
            {
              typedef typename boost::remove_all_extents<T>::type const result_type;
              result_type operator()( result_type &_t ) { return _t; }
            };
          template< class T1, class T2 >
            struct get_view< boost::fusion::joint_view<const T1, const T2> >
            {
              typedef boost::fusion::joint_view<const T1, const T2> const t_arg;
              typedef boost::fusion::single_view<t_arg const> result_type;
              result_type operator()( t_arg const &_t ) { return result_type(_t); }
            };
          template< class T >
            struct get_view< boost::fusion::single_view<const OptionalGroups<T> > >
            {
              typedef boost::fusion::single_view<const OptionalGroups<T> > t_arg;
              typedef typename OptionalGroups<T>::type result_type;
              result_type operator()( t_arg const &_t )
                { return boost::fusion::at_c<0>( boost::fusion::as_vector(_t) ).view; }
            };

        public:
          //! Resulting type.
          template< class T > struct result;
          //! Optional Group specialization.
          template< class T_THIS, class T1, class T2 >
            struct result<T_THIS(T1,T2)>
            {
              typedef OptionalGroups
                      <
                        typename join::result
                        <
                          join(typename get_view<T1>::result_type, 
                               typename get_view<T2>::result_type )
                        >::type
                      > type;
            };
          //! Optional Group specialization.
          template< class T_THIS, class T1, class T2 >
            struct result<T_THIS(T1 const&,T2 const&)>
            {
              typedef OptionalGroups
                      <
                        typename join::result
                        <
                          join(typename get_view<T1>::result_type, 
                               typename get_view<T2>::result_type )
                        >::type
                      > type;
            };
          //! Functor.
          template< class T1, class T2 >
            typename result<grouped_options(T1,T2)> :: type
              operator()( T1 const& _t1, T2 const& _t2 )
              {
                typedef typename result<grouped_options(T1,T2)>::type result_type;
                typename result_type::type result
                       ( 
                         join()
                         ( 
                           get_view<T1>()(_t1 ), 
                           get_view<T2>()(_t2 ) 
                         )
                       );
                result_type r(result);
                return r;
              }
      };


    } // namespace transform
  } // namespace load_n_save
} // namespace LaDa

namespace boost
{
  namespace fusion
  {
    namespace traits
    {
      template <class T>  struct is_view< LaDa::load_n_save::transform::OptionalGroups<T> >
      {
        typedef boost::mpl::false_ type;
      };
    }
  }
}
#endif

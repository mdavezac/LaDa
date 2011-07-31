#ifndef LADA_LNS_XPR_PARAMETERS_H
#define LADA_LNS_XPR_PARAMETERS_H

#include "LaDaConfig.h"

#include <boost/type_traits/is_same.hpp>
#include <boost/fusion/include/deref.hpp>
#include <boost/mpl/contains.hpp>
#include <boost/fusion/include/end.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/fusion/include/find_if.hpp>
#include <boost/mpl/placeholders.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/type_traits/remove_all_extents.hpp>
#include <boost/type_traits/remove_reference.hpp>
#include <boost/static_assert.hpp>


namespace LaDa 
{
  namespace load_n_save
  {
    namespace parameter
    {
      namespace details
      {
        struct no_type {};

        template<class T_TAG, class T_TYPE>
          struct assigned_parameter
          {
            typedef T_TAG t_Tag;
            typedef T_TYPE t_Type;
            t_Type type;
          };
      } 

      template< class T_TAG, class T_TYPE = details::no_type > struct parameter;

      template< class T_TAG, class T_TYPE >
        struct parameter<T_TAG, T_TYPE &>
        {
          typedef T_TAG t_Tag;
          details::assigned_parameter<T_TAG, T_TYPE&> operator=( T_TYPE &_type ) const
            { details::assigned_parameter<T_TAG, T_TYPE&> result = { _type }; return result; }
        };
      template< class T_TAG, class T_TYPE >
        struct parameter<T_TAG, T_TYPE const&>
        {
          typedef T_TAG t_Tag;
          details::assigned_parameter<T_TAG, T_TYPE const&> operator=( T_TYPE const &_type ) const
            { details::assigned_parameter<T_TAG, T_TYPE const&> result = { _type }; return result; }
        };
      template<class T_TAG>
        struct parameter<T_TAG, details::no_type>
        {
          typedef T_TAG t_Tag;
          template<class T_TYPE>
            details::assigned_parameter<T_TAG, T_TYPE const&> operator=( T_TYPE const& _type ) const
              { details::assigned_parameter<T_TAG, T_TYPE const&> result = { _type }; return result; }
          template<class T_TYPE>
            details::assigned_parameter<T_TAG, T_TYPE&> operator=( T_TYPE& _type ) const
              { details::assigned_parameter<T_TAG, T_TYPE&> result = { _type }; return result; }
        };
      template<class T_TAG>
        struct parameter<T_TAG, load_n_save::tags>
        {
          typedef T_TAG t_Tag;
          template<class T_TYPE>
            details::assigned_parameter<T_TAG, T_TYPE const&> operator=( T_TYPE const& _type ) const
              { details::assigned_parameter<T_TAG, T_TYPE const&> result = { _type }; return result; }
          template<class T_TYPE>
            details::assigned_parameter<T_TAG, T_TYPE&> operator=( T_TYPE& _type ) const
              { details::assigned_parameter<T_TAG, T_TYPE&> result = { _type }; return result; }
        };

      namespace details
      {
        template<class T0, class T1> bool got_param( T0 const& _t0, T1 &_t1, boost::mpl::true_ )
          { _t1 = _t0.type; return true; }
        template<class T0, class T1> bool got_param( T0 const& _t0, T1 &_t1, boost::mpl::false_ )
          { return false; }

        template<class T_TYPE, class T_DEFAULT >
          void get_param( T_TYPE &_type, T_DEFAULT const& _default, boost::mpl::false_ const& )
            { _type = _default; }

        template<class T_TYPE, class T_DEFAULT, class T_OTHER> 
          void get_param( T_TYPE &_type, T_DEFAULT const& _default, T_OTHER& _other )
            {  _type = _other; }

        template<class T_TAG>
          class GetParam
          {
            template< class T_FIRST, class T_LAST >
              struct is_last : public boost::fusion::result_of::equal_to<T_FIRST, T_LAST> {};
            template< class T_FIRST, class T_LAST, bool LAST > struct is_found;
            template< class T_FIRST, class T_LAST>
              struct is_found<T_FIRST, T_LAST, false> : public  boost::is_same
                   <
                     T_TAG,
                     typename boost::remove_reference
                     <
                       typename boost::fusion::result_of::deref<T_FIRST> :: type 
                     > :: type :: t_Tag
                   >  {};
            template< class T_FIRST, class T_LAST>
              struct is_found<T_FIRST, T_LAST, true> : public  boost::mpl::bool_<false> {};

            public:
              //! The resulting type.
              template<class T> struct result;
              //! The resulting type.
              template<class T_THIS, class T_FIRST, class T_LAST>
                struct result<T_THIS(T_FIRST, T_LAST)>
                {
                  //! Computes type depending on boolean.
                  template<class T, bool FOUND, bool LAST> struct what;
                  //! At last iterator.
                  template<class T, bool FOUND> struct what<T, FOUND, true>
                  {
                    //! Resulting type.
                    typedef boost::mpl::bool_<false> type;
                  };
                  //! Found tag.
                  template<class T>
                    struct what<T, true, false>  
                    {
                      //! Resulting type.
                      typedef typename boost::remove_reference
                              <
                                typename boost::fusion::result_of::deref
                                <
                                  T_FIRST
                                > :: type
                              > :: type :: t_Type type;
                    };
                  //! Found nothing, incrementing iterator.
                  template<class T>
                    struct what<T, false, false> 
                              : public result 
                                       <
                                         GetParam
                                         (
                                           typename boost::fusion::result_of::next<T_FIRST>::type,
                                           T_LAST
                                         )
                                        >  {};
                  //! The resulting type.
                  typedef typename what
                                   <
                                     void,
                                     is_found
                                     <
                                       T_FIRST, T_LAST,
                                       is_last<T_FIRST,T_LAST>::type::value
                                     >::type::value, 
                                     is_last<T_FIRST,T_LAST>::type::value
                                   > :: type type;
                };
  
              //! Functor.
              template<class T_FIRST, class T_LAST>
                typename result<GetParam(T_FIRST, T_LAST)>::type
                  operator()( T_FIRST const& _f, T_LAST const& _l ) const
                  {
                    return special_
                           ( 
                              _f, _l, 
                              typename is_found
                                <
                                  T_FIRST, T_LAST,
                                  is_last<T_FIRST, T_LAST>::type::value
                                >::type(),
                              typename is_last<T_FIRST, T_LAST>::type()  
                           ); 
                 }
             private:

              //! Is last iterator. 
              template<class T_FIRST, class T_LAST>
                typename result<GetParam(T_FIRST, T_LAST)>::type
                  special_( T_FIRST const&, T_LAST const&,
                            boost::mpl::false_ const&,
                            boost::mpl::true_ const&  ) const
                  { return boost::mpl::false_(); }
              //! Found tag.
              template<class T_FIRST, class T_LAST>
                typename result<GetParam(T_FIRST, T_LAST)>::type
                  special_( T_FIRST const& _f, T_LAST const&,
                            boost::mpl::true_ const&,
                            boost::mpl::false_ const&  ) const
                  {
                    //! Should not have more than one similar tag.
                    BOOST_STATIC_ASSERT
                    ((
                      boost::is_same
                      <
                        typename result
                        <
                          GetParam
                          (
                            typename boost::fusion::result_of::next<T_FIRST>::type,
                            T_LAST
                          )
                        > :: type, 
                        boost::mpl::bool_<false>
                      > :: type ::value
                    ));
                    return (*_f).type; 
                  }
              //! Found nothing, iterating.
              template<class T_FIRST, class T_LAST>
                typename result<GetParam(T_FIRST, T_LAST)>::type
                  special_( T_FIRST const& _f, T_LAST const& _l,
                            boost::mpl::false_ const&,
                            boost::mpl::false_ const&  ) const
                  { return operator()( boost::fusion::next(_f), _l ); }
          }; 
      }

      template<class T_TAG, class T_TYPE, class T_DEFAULT, class T_VECTOR >
        void get_param( T_TAG const& _t, T_TYPE &_type,
                        T_DEFAULT const& _default, T_VECTOR const& _vector )
        {
          typedef typename T_TAG :: t_Tag t_Tag;
          details::get_param
          (
            _type, _default, 
            details::GetParam<t_Tag>()
            (
              boost::fusion::begin(_vector),
              boost::fusion::end(_vector)
            )
          );
        }
      template<class T_TAG, class T_VECTOR >
        typename details::GetParam<typename T_TAG::t_Tag> :: template result
          <
            details::GetParam<typename T_TAG::t_Tag>
            (
              typename boost::fusion::result_of::begin<T_VECTOR>::type,
              typename boost::fusion::result_of::end<T_VECTOR>::type
            )
          > :: type get_param( T_TAG const&, T_VECTOR const& _vector )
          {
            return  details::GetParam<typename T_TAG::t_Tag>()
                    (
                      boost::fusion::begin(_vector),
                      boost::fusion::end(_vector)
                    );
          }

    }
  }
}

# endif

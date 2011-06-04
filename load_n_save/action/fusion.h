#ifndef BOOST_PP_IS_ITERATING
# ifndef LADA_LNS_XPR_FUSION_H
#   define LADA_LNS_XPR_FUSION_H
    
#   include "LaDaConfig.h"
    
#   include <boost/preprocessor/arithmetic/inc.hpp>
#   include <boost/preprocessor/repetition/enum_binary_params.hpp>
#   include <boost/preprocessor/repetition/enum_params.hpp>
    
#   include <boost/type_traits/remove_const.hpp>
#   include <boost/type_traits/remove_reference.hpp>
    
#   include<boost/tokenizer.hpp>
    
    
#   include <boost/fusion/algorithm/iteration/for_each.hpp>
#   include <boost/fusion/algorithm/iteration/accumulate.hpp>
#   include <boost/fusion/algorithm/transformation/transform.hpp>
#   include <boost/fusion/tuple.hpp>
    
    
#   include "../action/string_to_type.h"
#   include "../action/type_to_string.h"
#   include "../action/type_to_regex.h"
#   include "../action/type_to_regex.h"
#   include "../string_type.h"
    
    namespace LaDa
    {
      namespace load_n_save
      {
        namespace details
        {
          template<class T> struct remove_all :
            public boost::remove_const< typename boost::remove_reference<T>::type > {};

          template <typename First, typename Last>
            inline t_String type_to_regex(bool, boost::mpl::true_) { return ""; }
          template <typename First, typename Last>
            inline t_String type_to_regex(bool is_first, boost::mpl::false_) 
            {
              namespace bf = boost::fusion;
              typedef typename bf::result_of::next<First>::type t_next;
              typedef typename bf::result_of::equal_to<t_next, Last>::type t_test;
              typedef typename bf::result_of::deref<First>::type t_first;
              typedef typename remove_all<t_first>::type t_cvless;
              if(is_first) 
                return TypeToRegex<t_cvless>::apply() + type_to_regex<t_next, Last>(false, t_test());
              else  
                return "\\s+" + TypeToRegex<t_cvless>::apply()
                       + type_to_regex<t_next, Last>(false, t_test());
            }
          template <typename First, typename Last, typename T_TOKIT>
            inline bool string_to_type( First const &_first, Last const &_last, 
                                        T_TOKIT const &_tfirst, T_TOKIT const &_tlast,
                                        boost::mpl::true_)
              { return _tfirst == _tlast; }
          template <typename First, typename Last, typename T_TOKIT>
            inline bool string_to_type( First const &_first, Last const &_last, 
                                        T_TOKIT &_tfirst, T_TOKIT const &_tlast,
                                        boost::mpl::false_ ) 
            {
              if(_tfirst == _tlast) return false;
              namespace bf = boost::fusion;
              typedef typename bf::result_of::deref<First>::type t_first;
              typedef typename remove_all<t_first>::type t_cvless;
              if( not StringToType<t_cvless>::apply(*_tfirst, bf::deref(_first)) ) return false;
              typedef typename bf::result_of::next<First>::type t_next;
              typedef typename bf::result_of::equal_to<t_next, Last>::type t_test;
              ++_tfirst;
              return details::string_to_type(bf::next(_first), _last, _tfirst, _tlast, t_test()); 
            }

          template<class T> t_String type_to_string_(T const &_value)
            { return TypeToString<T>::apply(_value); }
          template <typename First, typename Last>
            inline t_String type_to_string( First const &_first, Last const &_last, 
                                            bool, boost::mpl::true_ )
              { return ""; }
          template <typename First, typename Last>
            inline t_String type_to_string( First const &_first, Last const &_last,
                                            bool is_first, boost::mpl::false_ )
            {
              namespace bf = boost::fusion;
              typedef typename bf::result_of::deref<First>::type t_first;
              typedef typename bf::result_of::next<First>::type t_next;
              typedef typename bf::result_of::equal_to<t_next, Last>::type t_test;
              typedef typename remove_all<t_first>::type t_cvless;
              if(is_first) 
                return   type_to_string_(bf::deref(_first))
                       + type_to_string(bf::next(_first), _last, false, t_test());
              else
                return   " " + type_to_string_(bf::deref(_first))
                       + type_to_string(bf::next(_first), _last, false, t_test());
            }
        }
      }
    }
#       define BOOST_PP_ITERATION_PARAMS_1 (3, (1, 4, <load_n_save/action/fusion.h>))
#       include BOOST_PP_ITERATE()
# endif
#else 
# define SIZE BOOST_PP_ITERATION()
  namespace LaDa
  {
    namespace load_n_save
    {
      //! Parses a fusion tuple.
      template<BOOST_PP_ENUM_PARAMS(SIZE, class T)>
        struct TypeToRegex<const boost::fusion::tuple<BOOST_PP_ENUM_PARAMS(SIZE, T)>, void>
        {
          //! Returns regex string.
          static t_String apply()
          { 
            namespace bf = boost::fusion;
            typedef boost::fusion::tuple<BOOST_PP_ENUM_PARAMS(SIZE, T)> t_type;
            typedef typename bf::result_of::begin<t_type>::type t_ifirst;
            typedef typename bf::result_of::end<t_type>::type t_ilast;
            typedef typename bf::result_of::equal_to<t_ifirst, t_ilast>::type t_test;
            return details::type_to_regex<t_ifirst, t_ilast>(true, t_test());
          }
        };
     
      //! Parses a fusion tuple.
      template<BOOST_PP_ENUM_PARAMS(SIZE, class T)>
        struct StringToType<const boost::fusion::tuple<BOOST_PP_ENUM_PARAMS(SIZE, T)>, void>
        {
          typedef const boost::fusion::tuple<BOOST_PP_ENUM_PARAMS(SIZE, T)> t_type;
          //! Functor.
          static bool apply( t_String const& _string, t_type &_value )
          {
            namespace bf = boost::fusion;
            boost::tokenizer<> tok(_string);
            boost::tokenizer<>::iterator i_first = tok.begin();
            boost::tokenizer<>::iterator i_last = tok.end();
            typedef typename bf::result_of::begin<t_type>::type t_ifirst;
            typedef typename bf::result_of::end<t_type>::type t_ilast;
            typedef typename bf::result_of::equal_to<t_ifirst, t_ilast>::type t_test;
            return details::string_to_type( bf::begin(_value), bf::end(_value),
                                            i_first, i_last, t_test() );
          }
        };

      template<BOOST_PP_ENUM_PARAMS(SIZE, class T)>
        struct TypeToString<const boost::fusion::tuple<BOOST_PP_ENUM_PARAMS(SIZE, T)>, void>
        {
          static t_String apply(boost::fusion::tuple<BOOST_PP_ENUM_PARAMS(SIZE, T)> const &_value)
          {
            namespace bf = boost::fusion;
            typedef boost::fusion::tuple<BOOST_PP_ENUM_PARAMS(SIZE, T)> t_type;
            typedef typename bf::result_of::begin<t_type>::type t_ifirst;
            typedef typename bf::result_of::end<t_type>::type t_ilast;
            typedef typename bf::result_of::equal_to<t_ifirst, t_ilast>::type t_test;
            return details::type_to_string(bf::begin(_value), bf::end(_value), true, t_test());
          }
        };
    }
  }
#endif

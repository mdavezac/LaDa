#ifndef BOOST_PP_IS_ITERATING
# ifndef   LADA_LNS_XPR_UTILITIES_H
#   define LADA_LNS_XPR_UTILITIES_H

#   include "LaDaConfig.h"

#   include <boost/preprocessor/iteration/iterate.hpp>
#   include <boost/preprocessor/repetition/enum_binary_params.hpp>
#   include <boost/preprocessor/repetition/enum_trailing_binary_params.hpp>
#   include <boost/preprocessor/facilities/intercept.hpp>
#   include <boost/fusion/include/vector.hpp>
#   include <boost/type_traits/is_const.hpp>
#   include <boost/type_traits/is_reference.hpp>
#   include <boost/type_traits/is_same.hpp>

#   include <opt/debug.h>
#   include "../string_type.h"
#   include "parameters.h"
#   include "section.h"
#   include "option.h"

    namespace LaDa 
    {
      namespace load_n_save
      {
        namespace details
        {
          struct help_tag {};
          struct tag_tag {};
          struct default_tag {};
          struct action_tag {};
        }
        const parameter::parameter<details::help_tag, std::string const&> help = {};
        const parameter::parameter<details::tag_tag, size_t const&> tag = {};
        const parameter::parameter<details::default_tag> default_ = {};
        const parameter::parameter<details::action_tag> action = {};
        inline xpr::Section section( t_String const& _name ) 
        {
          xpr::regular_data data;
          data.name = _name;
          data.help = "";
          data.tag = 0;
        
          xpr::Section section;
          section.set_data( data );
          return section;
        }
        
        inline xpr::Option option( t_String const& _name )
          { return xpr::Option(_name); }

        template< class T_TYPE > xpr::Section ext( T_TYPE &_a ) 
        {
          xpr::Section result;
          result.set_data( _a );
          return result;
        }

        template< class T_TYPE > xpr::Section ext( boost::shared_ptr<T_TYPE> const &_a ) 
        {
          xpr::Section result;
          result.set_data(_a);
          return result;
        }

        
#       define BOOST_PP_ITERATION_PARAMS_1 (3, (1, 4, <load_n_save/xpr/utilities.h>))
#       include BOOST_PP_ITERATE()
        
      }
    }
# endif
#else

# define SIZE BOOST_PP_ITERATION()
    template<BOOST_PP_ENUM_PARAMS(SIZE, class T)>
      xpr::Section section( t_String const& _name  
                            BOOST_PP_ENUM_TRAILING_BINARY_PARAMS(SIZE, T, const& _t) )
      {
        typedef boost::fusion::vector
           < BOOST_PP_ENUM_BINARY_PARAMS(SIZE, T, const& BOOST_PP_INTERCEPT) > t_Vector;
        t_Vector vec( BOOST_PP_ENUM_PARAMS(SIZE, _t) );
        xpr::regular_data data;
        data.name = _name;
        parameter::get_param( help, data.help, "", vec );
        parameter::get_param( tag, data.tag, 0, vec );

        xpr::Section section;
        section.set_data( data );
        return section;
      }

  template<BOOST_PP_ENUM_PARAMS(SIZE, class T)>
    xpr::Option option( t_String const& _name  
                        BOOST_PP_ENUM_TRAILING_BINARY_PARAMS(SIZE, T, const& _t) )
    {
      typedef boost::fusion::vector
         < BOOST_PP_ENUM_BINARY_PARAMS(SIZE, T, const& BOOST_PP_INTERCEPT) > t_Vector;
      t_Vector vec( BOOST_PP_ENUM_PARAMS(SIZE, _t) );
       
      // Static Assertions
      typedef typename parameter::details::GetParam<details::action_tag> :: template result
              <
                parameter::details::GetParam<details::action_tag>
                (
                  typename boost::fusion::result_of::begin<t_Vector>::type,
                  typename boost::fusion::result_of::end<t_Vector>::type
                )
              > :: type t_action;
      typedef typename parameter::details::GetParam<details::default_tag> :: template result
              <
                parameter::details::GetParam<details::default_tag>
                (
                  typename boost::fusion::result_of::begin<t_Vector>::type,
                  typename boost::fusion::result_of::end<t_Vector>::type
                )
              > :: type t_default;
      typedef typename boost::is_same<t_default, boost::mpl::bool_<false> > :: type no_default;
      typedef typename boost::is_same<t_action, boost::mpl::bool_<false> > :: type no_action;
      BOOST_STATIC_ASSERT(( no_default::value or (not no_action::value) ));
      BOOST_STATIC_ASSERT(( boost::is_reference<t_action>::type::value or no_action::value ));
      BOOST_STATIC_ASSERT(( boost::is_reference<t_default>::type::value or no_default::value ));
    // BOOST_STATIC_ASSERT
    // (( 
    //       (
    //         not boost::is_const
    //             < 
    //               typename boost::remove_reference<t_action>::type
    //             > :: type :: value
    //       )
    //    or no_action::value
    //    or action_::is_special_action
    //       < 
    //         typename boost::remove_const
    //         <
    //           typename boost::remove_reference<t_action> :: type
    //         > :: type
    //       > :: value
    // ));
      BOOST_STATIC_ASSERT
      (( 
         (
           boost::is_const
               < 
                 typename boost::remove_reference<t_default>::type
               > :: type :: value
         ) or no_default::value 
      ));


      xpr::Option op( _name );
      parameter::get_param( help, op.help, "", vec );
      parameter::get_param( tag, op.tag, 0, vec );
      op.set_action( parameter::get_param(action, vec), parameter::get_param(default_, vec) );

      // dynamic assertions.
      LADA_ASSERT( not (op.tag & load_n_save::required and (not no_default::value)),
                   "Required option " + op.name + " cannot have default value.\n" );
      LADA_ASSERT( not (op.tag & load_n_save::id and (not no_default::value)),
                   "Id option " + op.name + " cannot have default value.\n" );
      LADA_ASSERT( not (op.tag & load_n_save::id and op.tag & load_n_save::required),
                   "Option " + op.name + " cannot be both required and id.\n" );
      return op;
    }

#endif

#ifndef BOOST_PP_IS_ITERATING
# ifndef   LADA_LNS_XPR_UTILITIES_H
#   define LADA_LNS_XPR_UTILITIES_H

#   include "LaDaConfig.h"

#   include <sstream>

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
#   include "../tags.h"
#   include "../exceptions.h"
#   include "parameters.h"
#   include "section.h"
#   include "option.h"
#   include "../action/action.h"

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
          struct id_tag {};
        }
        const parameter::parameter<details::help_tag, std::string const&> help = {};
        const parameter::parameter<details::tag_tag, load_n_save::tags> tag = {};
        const parameter::parameter<details::default_tag> default_ = {};
        const parameter::parameter<details::action_tag> action = {};
        const parameter::parameter<details::id_tag> id = {};
        inline xpr::Section section( t_String const& _name ) 
        {
          xpr::regular_data data;
          data.name = _name;
          data.help = "";
          data.tag = required;
        
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
        xpr::Section ext(xpr::Section _a) { return _a; }

        template< class T_TYPE > xpr::Section ext( boost::shared_ptr<T_TYPE> const &_a ) 
        {
          xpr::Section result;
          result.set_data(_a);
          return result;
        }
        namespace details
        {
          void set_id_action(xpr::Option &_op, boost::mpl::bool_<false>) {}
          template<class T>
          void set_id_action(xpr::Option &_op, T const &_id)
          {
            _op.tag = idoption;
            std::ostringstream sstr; sstr << _id;
            _op.set_action(action_::IdAction(sstr.str()));
          }
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
        parameter::get_param( tag, data.tag, required, vec );
        if(data.tag == unavailable or data.tag == idoption)
          BOOST_THROW_EXCEPTION(error::invalid_section_tag());

        xpr::Section section;
        section.set_data( data );
        return section;
      }


  //! \brief Creates an xml option.
  //! \details This function takes a string as first argument: the name of the
  //!          xml option/attribute. A number of keyword arguments are allowed:
  //!            - tags=optional, if it is safe for the option not to be found
  //!                   in the xml input file.
  //!            - default_=value, where the value is assigned to the
  //!                   action if the option is not found in the inoput file.
  //!                   This keyword implies the one above.
  //!            - action=variable or special action, defines the action to be
  //!                   taken when the option is found in the document. If a
  //!                   variable, the value of the option is assigned to it.
  //!                   Note the this variable *cannot* be local. It must be
  //!                   defined when going out of scope. Special action allow
  //!                   for different possibilities, such as setting the
  //!                   numerical value of an enum type from strings, as in
  //!                   enum.h.
  //!            - help=string. Not yet supported. Reserved for future use.
  //!            - id=string. This allows XML tags to defined by their internal
  //!                   attributes. For instance all functionals could be found
  //!                   within separate <Functional> tags, but each ontaining a
  //!                   type attribute with a specific value, eg "vasp", "vff",
  //!                   etc. The parser will then check the existence of this
  //!                   attribute within the tag, and its value, before
  //!                   atempting to parse the xml it.
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
      typedef typename parameter::details::GetParam<details::id_tag> :: template result
              <
                parameter::details::GetParam<details::id_tag>
                (
                  typename boost::fusion::result_of::begin<t_Vector>::type,
                  typename boost::fusion::result_of::end<t_Vector>::type
                )
              > :: type t_id;
      typedef typename boost::is_same<t_default, boost::mpl::bool_<false> > :: type no_default;
      typedef typename boost::is_same<t_action, boost::mpl::bool_<false> > :: type no_action;
      typedef typename boost::is_same<t_id, boost::mpl::bool_<false> > :: type no_id;
      typedef typename boost::is_convertible<t_id, t_String > :: type has_id;
      BOOST_STATIC_ASSERT((    ((not no_id::value) and no_action::value and no_default::value)
                            or no_id::value ));
      BOOST_STATIC_ASSERT(( no_default::value or (not no_action::value) ));
      BOOST_STATIC_ASSERT(( boost::is_reference<t_action>::type::value or no_action::value ));
      BOOST_STATIC_ASSERT(( boost::is_reference<t_default>::type::value or no_default::value ));
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
      parameter::get_param( tag, op.tag, required, vec );
      if((not no_default::value) and op.tag & required ) op.tag = optional;
      if(op.tag == unavailable or op.tag == idoption)
        BOOST_THROW_EXCEPTION(error::invalid_section_tag());
      if((not no_id::value))
        details::set_id_action(op, parameter::get_param(id, vec));
      else
        op.set_action( parameter::get_param(action, vec), parameter::get_param(default_, vec) );

      // dynamic assertions.
      LADA_ASSERT( not (op.tag & load_n_save::required and (not no_default::value)),
                   "Required option " + op.name + " cannot have default value.\n" );
      LADA_ASSERT( not (op.tag & load_n_save::idoption and (not no_default::value)),
                   "Id option " + op.name + " cannot have default value.\n" );
      LADA_ASSERT( not (op.tag & load_n_save::idoption and op.tag & load_n_save::required),
                   "Option " + op.name + " cannot be both required and id.\n" );
      return op;
    }

#endif

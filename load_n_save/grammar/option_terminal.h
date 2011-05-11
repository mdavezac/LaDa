//
//  Version: $Id: terminals.h 1200 2009-06-22 05:18:21Z davezac $
//

#if !BOOST_PP_IS_ITERATING
# ifndef _LADA_LOADNSAVE_GRAMMAR_OPTION_TERMINAL_H_
#   define _LADA_LOADNSAVE_GRAMMAR_OPTION_TERMINAL_H_
    
#   ifdef HAVE_CONFIG_H
#   include <config.h>
#   endif
    
#   ifdef _LADADEBUG
#    include <iostream>
#   endif
    
#   include "terminals.h"
#   include "tags.h"
    
    namespace LaDa 
    {
      namespace load_n_save
      {
        namespace grammar
        {
          namespace proto = boost::proto;
    
          namespace details
          {
            struct option_function_tag {};
            struct name_tag {};
            struct tag_tag {};
            struct help_tag {};
            struct default_tag {};
            struct action_tag {};
          }
    
          proto::terminal<details::tag_tag>::type const tags = {{}};
          proto::terminal<details::help_tag>::type const help = {{}};
          proto::terminal<details::default_tag>::type const default_ = {{}};
          proto::terminal<details::action_tag>::type const action = {{}};
    
          struct tag_grammar : proto::assign
            < 
              proto::terminal<details::tag_tag>,
              proto::or_
              <
                proto::terminal<tags::option::option>,
                proto::terminal<tags::section::section>,
                proto::terminal<int>
              >
            > {};
    
          struct help_grammar : proto::assign
            < 
              proto::terminal<details::help_tag>,
              string_terminal
            > {};
    
          struct default_grammar : proto::assign
            < 
              proto::terminal<details::default_tag>,
              proto::terminal<proto::_ const &>
            > {};
    
          struct action_grammar : proto::assign
            < 
              proto::terminal<details::action_tag>,
              proto::or_
              <
                proto::terminal<proto::_ const &>,
                proto::terminal<proto::_&>
              >
            > {};
    
          struct option_function : proto::function
            <
              proto::terminal<details::option_function_tag>::type,
              string_terminal,
              proto::vararg
              <
                proto::or_
                <
                  tag_grammar,
                  help_grammar,
                  default_grammar,
                  action_grammar
                >
              >
            > {};
    
        } // namespace grammar
    
      } // namespace load_n_save
    
    } // namespace LaDa

#   define BOOST_PP_ITERATION_PARAMS_1 (3, (1, 5, <load_n_save/grammar/option_terminal.h>))
#   include BOOST_PP_ITERATE()
# endif

#else

# define SIZE BOOST_PP_ITERATION()
  namespace LaDa 
  {
    namespace load_n_save
    {
      namespace grammar
      {

#       ifdef LADA_ENUM_CREF_PARAMS
#         error LADA_ENUM_CREF_PARAMS macro already exists.
#       endif
#       ifdef LADA_ENUM_CREF_BINARY_PARAMS
#         error LADA_ENUM_CREF_BINARY_PARAMS macro already exists.
#       endif
#       ifdef LADA_ENUM_ARGS
#         error LADA_ENUM_ARGS macro already exists.
#       endif
#       define LADA_ENUM_CREF_PARAMS(z, N, data) T ## N const&
#       define LADA_ENUM_CREF_BINARY_PARAMS(z, N, data) T ## N const& _t ## N
#       define LADA_ENUM_ARGS(z, N, data) boost::cref( _t ## N )
          
        template< BOOST_PP_ENUM_PARAMS(SIZE, class T) >
          typename proto::result_of::make_expr
          <
            proto::tag::function,
            Domain,
            details::option_function_tag const,
            BOOST_PP_ENUM(SIZE, LADA_ENUM_CREF_PARAMS, )
          > :: type option( BOOST_PP_ENUM(SIZE, LADA_ENUM_CREF_BINARY_PARAMS, ) )
          {
            return proto::make_expr<proto::tag::function, Domain>
                   ( 
                     details::option_function_tag()
                     BOOST_PP_ENUM_TRAILING( SIZE, LADA_ENUM_ARGS, )
                   );
          }
        
#       undef LADA_ENUM_CREF_PARAMS
#       undef LADA_ENUM_CREF_BINARY_PARAMS
#       undef LADA_ENUM_ARGS


      } // namespace grammar
    } // namespace load_n_save
  } // namespace LaDa


#endif

//
//  Version: $Id: section_terminal.h 1211 2009-07-04 02:03:24Z davezac $
//

#if !BOOST_PP_IS_ITERATING
# ifndef _LADA_LOADNSAVE_GRAMMAR_SECTION_TERMINAL_H_
#   define _LADA_LOADNSAVE_GRAMMAR_SECTION_TERMINAL_H_
    
#   ifdef HAVE_CONFIG_H
#   include <config.h>
#   endif
    
#   include <boost/preprocessor/iteration/iterate.hpp>
#   include <boost/preprocessor/repetition/enum.hpp>
#   include <boost/preprocessor/repetition/enum_params.hpp>
    
#   include "terminals.h"
#   include "option_terminal.h"
    
    namespace LaDa 
    {
      namespace load_n_save
      {
        namespace grammar
        {
          namespace proto = boost::proto;
    
          namespace details
          {
            struct section_function_tag {};
          }
    
    
          struct section_function : proto::function
            <
              proto::terminal<details::section_function_tag>,
              string_terminal,
              proto::vararg
              <
                proto::or_
                <
                  tag_grammar,
                  help_grammar
                >
              >
            > {};
        } // namespace grammar
      } // namespace load_n_save
    } // namespace LaDa

#   define BOOST_PP_ITERATION_PARAMS_1 (3, (1, 3, <load_n_save/grammar/section_terminal.h>))
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
            details::section_function_tag const,
            BOOST_PP_ENUM(SIZE, LADA_ENUM_CREF_PARAMS, )
          > :: type section( BOOST_PP_ENUM(SIZE, LADA_ENUM_CREF_BINARY_PARAMS, ) )
          {
            return proto::make_expr<proto::tag::function, Domain>
                   ( 
                     details::section_function_tag()
                     BOOST_PP_ENUM_TRAILING( SIZE, LADA_ENUM_ARGS, )
                   );
          }
        
#       undef LADA_ENUM_CREF_PARAMS
#       undef LADA_ENUM_CREF_BINARY_PARAMS
#       undef LADA_ENUM_ARGS


      } // namespace grammar
    } // namespace load_n_save
  } // namespace LaDa


# endif

//
//  Version: $Id: terminals.h 1205 2009-06-24 00:07:30Z davezac $
//

#ifndef _LADA_LOADNSAVE_GRAMMAR_TERMINALS_H_
#define _LADA_LOADNSAVE_GRAMMAR_TERMINALS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef _LADADEBUG
# include <iostream>
#endif
#include <boost/fusion/include/at_c.hpp>

#include <boost/proto/tags.hpp>
#include <boost/proto/traits.hpp>
#include <boost/proto/literal.hpp>

#include <opt/types.h>

#include "dsel.h"
#include "tags.h"
#include "external_type.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace grammar
    {
      namespace proto = boost::proto;


      namespace details
      {
#       ifdef LADA_BODY
#         error LADA_BODY already defined.
#       endif
#       ifdef LADA_PRINT
#         error LADA_PRINT already defined.
#       endif
       
#       define LADA_BODY \
          char const *value; \
          size_t const tag; \
          char const *help; \
          operator char const*() const { return value; } \
      
#       define LADA_PRINT( what )  \
          inline std::ostream & operator<<( std::ostream& _stream, const what & _i ) \
              { return _stream << _i.value; }
        //! An option data structure.
        struct option_terminal { LADA_BODY }; 
        //! A section data structure. 
        struct section_terminal { LADA_BODY };

        LADA_PRINT( option_terminal )
        LADA_PRINT( section_terminal )

#       undef LADA_BODY
#       undef LADA_PRINT

        //! Tag for dynamic expression terminals.
        struct dynamic_expression {};
      }

      //! Proto option terminal.
      typedef proto::terminal<details::option_terminal> option_terminal;
      //! Proto option terminal.
      typedef proto::terminal<details::section_terminal> section_terminal;
      //! Dynamic expression terminals.
      typedef proto::terminal<details::dynamic_expression> dynamic_expression;

     
      //! A terminal to match a string.
      struct string_terminal : public proto::or_
        <
          proto::terminal< const char* >,
          proto::terminal< char const[proto::N] >,
          proto::terminal< std::string const >
        > {};
     
      //! A terminal to hold a (non-constant) reference.
      template< class  T > Dsel<typename proto::terminal<T&>::type> const lit(T &t)
        {
          Dsel< typename proto::terminal<T &>::type > op = {t};
          return op;
        }
      //! A terminal to hold a nullary constant callable.
      template< class  T > 
        Dsel<typename proto::terminal< T const& >::type> const lit(const T &t)
        {
          Dsel< typename proto::terminal<T const&>::type > op = {t};
          return op;
        }
     
//     //! A terminal to hold an option.
//     inline const Dsel< option_terminal::type >
//       option( const char *t, size_t _tag = tags::section::default_ )
//       {
//         Dsel< option_terminal::type > op = {{t, _tag}};
//         return op;
//       }
//     //! A terminal to hold an option.
//     inline const Dsel< section_terminal::type >
//       section( const char *t, size_t _tag = tags::option::default_ )
//       {
//         Dsel< section_terminal::type > sec = {{t, _tag}};
//         return sec;
//       }

    } // namespace grammar

  } // namespace load_n_save

} // namespace LaDa


#endif

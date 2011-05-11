//
//  Version: $Id: values.h 1205 2009-06-24 00:07:30Z davezac $
//

#ifndef _LADA_LOADNSAVE_GRAMMAR_VALUES_H_
#define _LADA_LOADNSAVE_GRAMMAR_VALUES_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef _LADADEBUG
# include <iostream>
#endif
#include <boost/fusion/include/at_c.hpp>

#include <boost/proto/expr.hpp>
#include <boost/proto/extends.hpp>
#include <boost/proto/tags.hpp>
#include <boost/proto/traits.hpp>
#include <boost/proto/operators.hpp>
#include <boost/proto/domain.hpp>


#include <opt/types.h>

#include "terminals.h"
#include "option_terminal.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace grammar
    {
      namespace proto = boost::proto;
      
      template< class T_TYPE >
        struct TypedValue : public proto::or_
          <
            proto::terminal< T_TYPE >,
            proto::equal_to< string_terminal, proto::terminal< T_TYPE > >
          > {};

      template< class T_TYPE >
        struct TypedValues : public proto::or_
          <
            TypedValue<T_TYPE>,
            proto::logical_or< TypedValues<T_TYPE>, TypedValues<T_TYPE> >
          > {};

      //! All possible assignements for an option.
      struct Value : public  TypedValue< proto::_ > {};

      //! All possible assignements for an option.
      struct Values : public  TypedValues< proto::_ > {};
      
    } // namespace grammar

  } // namespace load_n_save

} // namespace LaDa


#endif

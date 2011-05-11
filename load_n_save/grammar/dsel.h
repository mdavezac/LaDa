//
//  Version: $Id: dsel.h 1231 2009-07-17 05:12:39Z davezac $
//

#ifndef _LADA_LOADNSAVE_DSEL_H_
#define _LADA_LOADNSAVE_DSEL_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/proto/expr.hpp>
#include <boost/proto/extends.hpp>

#include <opt/types.h>
#include "external_type.h"
#include "action_type.h"


namespace LaDa 
{
  namespace load_n_save
  {
    namespace grammar
    {

      namespace proto = boost::proto;
  
      //! \cond 
      //  The domain for the load_n_save dsel.
      struct Domain;
      //! \cond
  
      //! Domain specific embeded language. Pure grammar.
      template< typename T_EXPR >
        struct Dsel
        {
            BOOST_PROTO_BASIC_EXTENDS(T_EXPR, Dsel<T_EXPR>, Domain)
            BOOST_PROTO_EXTENDS_ASSIGN()
//           BOOST_PROTO_EXTENDS_SUBSCRIPT()
        };
      // includes specializations from other files.
#     ifdef FROM_LADA_DSEL
#       error FROM_LADA_DSEL macro already defined.
#     endif
#     define FROM_LADA_DSEL

#       include "external_type.h"
#       include "action_type.h"

#     undef FROM_LADA_DSEL
    } // namespace grammar

  } // namespace load_n_save

} // namespace LaDa


#endif

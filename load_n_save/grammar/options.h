//
//  Version: $Id: options.h 1215 2009-07-08 00:28:08Z davezac $
//

#ifndef _LADA_LOADNSAVE_GRAMMAR_OPTIONS_H_
#define _LADA_LOADNSAVE_GRAMMAR_OPTIONS_H_

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
#include "values.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace grammar
    {
      namespace proto = boost::proto;
      
      //! Eventually, and (actioned) option can be ornated in the following flavors.
      struct OrnateOption : public proto::or_
        <
          option_function,
          proto::assign< option_function, TypedValues<proto::_> >
        > {};
      
      //! Options are linked together through comma.
      struct OrnateOptions : public proto::or_
        <
          OrnateOption,
          proto::plus< OrnateOptions, OrnateOptions >,
          proto::logical_or< OrnateOptions, OrnateOptions >
        > {};

    } // namespace grammar
 

  } // namespace load_n_save

} // namespace LaDa


#endif

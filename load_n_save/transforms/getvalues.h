//
//  Version: $Id: getvalues.h 1211 2009-07-04 02:03:24Z davezac $
//

#ifndef _LADA_LOADNSAVE_TRANSFORMS_GETVALUES_H_
#define _LADA_LOADNSAVE_TRANSFORMS_GETVALUES_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/proto/core.hpp>
#include <boost/proto/transform.hpp>
#include <boost/proto/transform/when.hpp>
#include <boost/proto/matches.hpp>

#include "../grammar/options.h"
#include "protect.h"
#include "join.h"
#include "false.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace transform
    {
      namespace proto = boost::proto;

      struct GetValues : public proto::or_
        <
          proto::when< grammar::option_function, false_() >,
          proto::when< proto::terminal< proto::_ >, protect( proto::_value ) >,
          proto::when
          <
            proto::equal_to< grammar::string_terminal, proto::terminal< proto::_ > >,
            protect( join( GetValues(proto::_left), GetValues(proto::_right) ) )
          >,
          proto::when
          < 
            proto::logical_or< grammar::Values, grammar::Values >,
            join( GetValues(proto::_left), GetValues(proto::_right) )
          >,
          proto::when
          <
            proto::assign< grammar::option_function, grammar::Values >,
            GetValues(proto::_right)
          >,
          proto::otherwise< false_() >
        > {};

    } // namespace transform
 
  } // namespace load_n_save

} // namespace LaDa


#endif

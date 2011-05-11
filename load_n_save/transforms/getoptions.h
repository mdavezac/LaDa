//
//  Version: $Id: getoptions.h 1226 2009-07-13 06:28:01Z davezac $
//

#ifndef _LADA_LOADNSAVE_TRANSFORMS_SECTIONS_H_
#define _LADA_LOADNSAVE_TRANSFORMS_SECTIONS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/proto/core.hpp>
#include <boost/proto/transform.hpp>
#include <boost/proto/transform/when.hpp>
#include <boost/proto/matches.hpp>

#include "../grammar/options.h"
#include "group_optionals.h"
#include "getcontent.h"
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


//     namespace details
//     {
//       struct GetOptionalSequence;
//     }

      struct GetOptions : public proto::or_
        <
          proto::when< grammar::OrnateOption, protect( proto::_ ) >,
          proto::when
          <
            proto::plus< grammar::OrnateOptions, grammar::OrnateOption >,
            join( GetOptions(proto::_left), protect(proto::_right) )
          >,
          proto::when
          <
            proto::plus< grammar::OrnateOptions, grammar::OrnateOptions >,
            join( GetOptions(proto::_left), GetOptions(proto::_right) )
          >,
          proto::when
          <
            proto::plus< grammar::Content, grammar::OrnateSection >,
            GetOptions(proto::_left) 
          >,
          proto::when
          <
            proto::logical_or< grammar::OrnateOptions, grammar::OrnateOptions >,
            protect
            ( 
              transform::grouped_options
              ( 
                GetOptions(proto::_left), 
                GetOptions(proto::_right)
              )
            )
          >,
          proto::otherwise< false_() >
        > {};

//     namespace details
//     {
//       struct GetOptionalSequence : public proto::or_
//              <
//                proto::when
//                <
//                  proto::logical_or< grammar::OrnateOptions, grammar::OrnateOptions>,
//                  join( protect(GetOptions(proto::_left)), join(GetOptions(proto::_right_)) )
//                >
//              >  {};
//     }

    } // namespace transform
 
  } // namespace load_n_save

} // namespace LaDa


#endif

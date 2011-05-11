//
//  Version: $Id: getsections.h 1226 2009-07-13 06:28:01Z davezac $
//

#ifndef _LADA_LOADNSAVE_TRANSFORMS_GETSECTIONS_H_
#define _LADA_LOADNSAVE_TRANSFORMS_GETSECTIONS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/proto/core.hpp>
#include <boost/proto/transform.hpp>
#include <boost/proto/transform/when.hpp>

#include "../grammar/sections.h"
#include "protect.h"
#include "join.h"
#include "options.h"
#include "false.h"
#include "group_optionals.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace transform
    {
      namespace proto = boost::proto;

      struct GetSections : public proto::or_
        <
          proto::when< grammar::OrnateSection, protect( proto::_ ) >,
          proto::when
          <
            proto::plus< grammar::OrnateOptions, grammar::OrnateSection >,
            GetSections(proto::_right) 
          >,
          proto::when
          <
            proto::plus< grammar::Content, grammar::OrnateSection >,
            join( GetSections(proto::_left), protect(proto::_right)  )
          >,
          proto::when
          <
            proto::plus< grammar::OrnateSections, grammar::OrnateSection >,
            protect( grouped_options( GetSections(proto::_left), protect(proto::_right) ) )
          >,
          proto::when
          <
            proto::plus< grammar::OrnateSections, grammar::OrnateSections >,
            join( GetSections(proto::_left), GetSections(proto::_right) )
          >,
          proto::when
          <
            proto::logical_or< grammar::OrnateSections, grammar::OrnateSections >,
            protect
            ( 
              transform::grouped_options
              ( 
                GetSections(proto::_left), 
                GetSections(proto::_right)
              )
            )
          >,
          proto::otherwise< false_() >
        > {};

  
    } // namespace transform
 
  } // namespace load_n_save

} // namespace LaDa


#endif

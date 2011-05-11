//
//  Version: $Id: getcontent.h 1170 2009-06-09 04:17:08Z davezac $
//

#ifndef _LADA_LOADNSAVE_TRANSFORMS_GETCONTENT_H_
#define _LADA_LOADNSAVE_TRANSFORMS_GETCONTENT_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/proto/core.hpp>
#include <boost/proto/transform.hpp>
#include <boost/proto/transform/when.hpp>

#include "../grammar/sections.h"
#include "join.h"
#include "false.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace transform
    {
      namespace proto = boost::proto;

      //! Returns content of a section (e.g options + subsections )
      struct GetContent : public proto::or_
        <
          proto::when
          <
            proto::subscript< grammar::OrnateSection, proto::terminal<proto::_ const&> >,
            GetContent( proto::_left )
          >,
          proto::when
          <
            proto::assign< grammar::OrnateSection, grammar::Content >,
            proto::_right
          >,
          proto::otherwise< boost::mpl::bool_<false>() >
        > {};
    } // namespace transform
 
  } // namespace load_n_save

} // namespace LaDa


#endif

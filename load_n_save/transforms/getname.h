//
//  Version: $Id: getname.h 1211 2009-07-04 02:03:24Z davezac $
//

#ifndef _LADA_LOADNSAVE_TRANSFORMS_GETNAME_H_
#define _LADA_LOADNSAVE_TRANSFORMS_GETNAME_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/proto/core.hpp>
#include <boost/proto/transform.hpp>
#include <boost/proto/transform/when.hpp>

#include "../grammar/sections.h"
#include "../grammar/options.h"
#include "terminals.h"
#include "false.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace transform
    {
      namespace proto = boost::proto;
      namespace fusion = boost::fusion;

      //! Returns name of an option or a section.
      struct GetName : public  proto::or_
        <
          proto::when< grammar::option_function, proto::_value(proto::_child_c<1>) >,
          proto::when< grammar::section_function, proto::_value(proto::_child_c<1>) >,
          proto::when
          <
            proto::assign< grammar::section_function, grammar::Content >,
            GetName( proto::_left )
          >,
          proto::when
          <
            proto::assign< grammar::option_function, grammar::Values >,
            GetName( proto::_left )
          >,
          proto::otherwise< boost::mpl::bool_<false>() >
        > {};

    } // namespace transform
 
  } // namespace load_n_save

} // namespace LaDa


#endif

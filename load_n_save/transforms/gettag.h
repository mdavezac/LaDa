//
//  Version: $Id: gettag.h 1211 2009-07-04 02:03:24Z davezac $
//

#ifndef _LADA_LOADNSAVE_TRANSFORMS_GETTAG_H_
#define _LADA_LOADNSAVE_TRANSFORMS_GETTAG_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/proto/core.hpp>
#include <boost/proto/transform.hpp>
#include <boost/proto/transform/when.hpp>

#include "../grammar/sections.h"
#include "../grammar/options.h"
#include "false.h"
#include "find_named_argument.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace transform
    {
      namespace proto = boost::proto;

      //! Returns tags of section or option.
      struct GetTag : public proto::or_ 
        < 
          proto::when 
          <  
            grammar::option_function,  
            details::ReturnValue< 2, grammar::tag_grammar, boost::mpl::int_<0>  >()  
          >, 
          proto::when 
          <  
            grammar::section_function,  
            details::ReturnValue< 2, grammar::tag_grammar, boost::mpl::int_<0>  >()  
          >, 
          proto::when 
          < 
            proto::assign< grammar::section_function, grammar::Content >, 
            GetTag( proto::_left ) 
          >, 
          proto::when 
          < 
            proto::assign< grammar::option_function, grammar::Values >, 
            GetTag( proto::_left ) 
          >, 
          proto::otherwise< boost::mpl::int_<0> () > 
        > {};

    } // namespace transform
 
  } // namespace load_n_save

} // namespace LaDa


#endif

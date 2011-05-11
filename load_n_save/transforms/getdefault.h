//
//  Version: $Id: getdefault.h 1212 2009-07-04 04:36:19Z davezac $
//

#ifndef _LADA_LOADNSAVE_TRANSFORMS_GETDEFAULT_H_
#define _LADA_LOADNSAVE_TRANSFORMS_GETDEFAULT_H_

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
      struct GetDefault : public proto::or_ 
        < 
          proto::when 
          <  
            grammar::option_function,  
            details::ReturnValue< 2, grammar::default_grammar, boost::mpl::bool_<false>  >()  
          >, 
          proto::when 
          < 
            proto::assign< grammar::option_function, grammar::Values >, 
            GetDefault( proto::_left ) 
          >, 
          proto::otherwise< boost::mpl::bool_<false>() > 
        > {};

    } // namespace transform
 
  } // namespace load_n_save

} // namespace LaDa


#endif

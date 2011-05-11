//
//  Version: $Id: gethelp.h 1211 2009-07-04 02:03:24Z davezac $
//

#ifndef _LADA_LOADNSAVE_TRANSFORMS_GETHELP_H_
#define _LADA_LOADNSAVE_TRANSFORMS_GETHELP_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "find_named_argument.h"
#include "false.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace transform
    {
      namespace proto = boost::proto;
      namespace details
      {
        struct empty_string
        {
          operator char const*() { return ""; }
        };
      };

      //! Returns help string of section or option.
      struct GetHelp : public proto::or_
        < 
          proto::when 
          <  
            grammar::option_function,  
            details::ReturnValue< 2, grammar::help_grammar, details::empty_string >()  
          >, 
          proto::when 
          <  
            grammar::section_function,  
            details::ReturnValue< 2, grammar::help_grammar, details::empty_string >()  
          >, 
          proto::when 
          < 
            proto::assign< grammar::section_function, grammar::Content >, 
            GetHelp( proto::_left ) 
          >, 
          proto::when 
          < 
            proto::assign< grammar::option_function, grammar::Values >, 
            GetHelp( proto::_left ) 
          >, 
          proto::otherwise< details::empty_string() > 
        > {};
  
    } // namespace transform
 
  } // namespace load_n_save

} // namespace LaDa


#endif

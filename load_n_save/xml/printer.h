//
//  Version: $Id: printer.h 1167 2009-06-07 23:22:42Z davezac $
//

#ifndef _LADA_LOADNSAVE_XML_FORMAT_H_
#define _LADA_LOADNSAVE_XML_FORMAT_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "../tree/tree.h"

namespace LaDa 
{
  namespace load_n_save
  { 
    namespace xml
    {
      //! Prints a parsed tree in XML format.
      std::ostream& print( std::ostream &_stream, tree::Base const& _base,
                           std::string const &_indent = "", std::string const &_char = "  " );
      //! Prints a parsed tree in XML format.
      std::ostream& print( std::ostream &_stream, tree::Section const& _sec,
                           std::string const &_indent = "", std::string const &_char = "  " );
    } // namespace xml
  } // namespace load_n_save
} // namespace LaDa

#endif

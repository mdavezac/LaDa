//
//  Version: $Id: parser.h 1226 2009-07-13 06:28:01Z davezac $
//

#ifndef _LADA_LOADNSAVE_XML_PRINTER_H_
#define _LADA_LOADNSAVE_XML_PRINTER_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "../tree/section.h"

namespace LaDa 
{
  namespace load_n_save
  { 
    namespace xml
    {
//     //! Returns a section containing parsed file.
//     boost::shared_ptr<tree::Base> parse( std::string const &_path );
      //! Returns a section containing parsed from a string.
      boost::shared_ptr<tree::Base> parse( std::string const &_text );

    } // namespace xml
  } // namespace load_n_save
} // namespace LaDa

#endif

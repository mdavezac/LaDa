#ifndef LADA_LOADNSAVE_XML_PRINTER_H
#define LADA_LOADNSAVE_XML_PRINTER_H

#include "LaDaConfig.h"

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

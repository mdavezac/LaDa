#ifndef LADA_LOADNSAVE_XML_FORMAT_H
#define LADA_LOADNSAVE_XML_FORMAT_H

#include "LaDaConfig.h"
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

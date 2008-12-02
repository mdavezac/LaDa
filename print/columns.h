//
//  Version: $Id: base.h 844 2008-11-08 01:22:54Z davezac $
//
#ifndef _LADA_PRINT_COLUMNS_H_
#define _LADA_PRINT_COLUMNS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>

namespace LaDa
{
  namespace Print
  {
     //! Two column output. Taken from boost::program_options
     void format_one( std::ostream& _stream, const std::string& _first, 
                      const std::string& _second,
                      size_t _first_column_width, size_t _line_length);
  } // namespace Print
} // namespace LaDa
#endif 

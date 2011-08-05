#ifndef LADA_LNS_EXCEPTIONS_H
#define LADA_LNS_EXCEPTIONS_H

#include "LaDaConfig.h"

#include <iostream>

#include <boost/throw_exception.hpp>
#include <boost/exception/error_info.hpp>
#include <boost/exception/info.hpp>


namespace LaDa
{
  namespace load_n_save
  {
    namespace error
    {
      //! Internal error.
      struct internal: virtual boost::exception, virtual std::exception {};

      //! Thrown when a section is explicitely given an invalid tag.
      struct invalid_section_tag : virtual internal {};
      //! Thrown when an option is explicitely given an invalid tag.
      struct invalid_option_tag : virtual internal {};

      //! Root of section exceptions.
      struct section: virtual boost::exception, virtual std::exception {};
      //! Thrown when a required section is not found.
      struct section_not_found: virtual section {};
      //! Thrown when a section is not found.
      struct required_section_not_found: virtual section_not_found {}; 
      //! Thrown when to few of a given sections are found, eg <Atom/>.
      struct too_few_sections: virtual section_not_found {}; 
      //! Thrown when parsing did not succeed.
      struct section_parse_error: virtual section {};

      //! Cannot merge a recurrent expresion, such as push_back.
      struct cannot_merge_recurrent : virtual internal {};
      //! Cannot load or save an empty shared ptr. 
      struct empty_shared_ptr : virtual internal {};

      //! Root of option exceptions.
      struct option: virtual boost::exception, virtual std::exception {};
      //! Thrown when a required option is not found.
      struct option_not_found: virtual option {};
      //! Thrown when a option is not found.
      struct required_option_not_found: virtual option_not_found {}; 
      //! Thrown when parsing did not succeed.
      struct option_parse_error: virtual option {};
      //! Thrown when parsing did not succeed.
      struct enum_transcript_error: virtual option_parse_error {};

      //! Root of excepetions when parsing xml strings.
      struct parser_error : virtual boost::exception, virtual std::exception {};
      //! Thrown when end of tag cannot be found.
      struct no_end_tag : virtual parser_error {};
      //! Thrown when option cannot be parsed.
      struct unparsable_option : virtual parser_error {};
      //! Thrown when no xml tag is found beyond a given line.
      struct no_xml_tag : virtual parser_error {};
      //! Thrown when input cannot be parsed.
      struct unparsable_text : virtual parser_error {};

      //! Name of the section
      typedef boost::error_info<struct lns_sec_name,std::string> section_name;
      //! Name of the section
      typedef boost::error_info<struct lns_op_name,std::string> option_name;
      //! Line number.
      typedef boost::error_info<struct lns_line_number, size_t> line_number;
    }
  }
}

#endif

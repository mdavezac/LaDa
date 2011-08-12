#ifndef LADA_LNS_EXCEPTIONS_H
#define LADA_LNS_EXCEPTIONS_H

#include "LaDaConfig.h"

#include <root_exceptions.h>


namespace LaDa
{
  namespace error
  {
    //! Thrown when a section is explicitely given an invalid tag.
    struct invalid_section_tag : virtual internal {};
    //! Thrown when an option is explicitely given an invalid tag.
    struct invalid_option_tag : virtual internal {};

    //! Root of section exceptions.
    struct section: virtual root {};
    //! Thrown when a required section is not found.
    struct section_not_found: virtual section, virtual input {};
    //! Thrown when a section is not found.
    struct required_section_not_found: virtual section_not_found {}; 
    //! Thrown when to few of a given sections are found, eg <Atom/>.
    struct too_few_sections: virtual section_not_found {}; 
    //! Thrown when parsing did not succeed.
    struct section_parse_error: virtual section, virtual input {};

    //! Cannot merge a recurrent expresion, such as push_back.
    struct cannot_merge_recurrent : virtual internal {};
    //! Cannot load or save an empty shared ptr. 
    struct empty_shared_ptr : virtual internal {};

    //! Root of option exceptions.
    struct option: virtual root {};
    //! \brief Thrown when a option is not found.
    //! \details This is not an error, merely a way to unwind the stack.
    struct option_not_found: virtual option {};
    //! Thrown when a required option is not found.
    struct required_option_not_found: virtual option_not_found, virtual input {}; 
    //! Thrown when parsing did not succeed.
    struct option_parse_error: virtual option, virtual input {};
    //! Thrown when parsing did not succeed.
    struct enum_transcript_error: virtual option_parse_error {};

    //! Root of excepetions when parsing xml strings.
    struct parser_error : virtual input {};
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

#endif

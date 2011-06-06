#ifndef LADA_LNS_EXCEPTIONS_H_
#define LADA_LNS_EXCEPTIONS_H_

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

      //! Root of option exceptions.
      struct option: virtual boost::exception, virtual std::exception {};
      //! Thrown when a required option is not found.
      struct option_not_found: virtual option {};
      //! Thrown when a option is not found.
      struct required_option_not_found: virtual option_not_found {}; 
      //! Thrown when parsing did not succeed.
      struct option_parse_error: virtual option {};

      //! Name of the section
      typedef boost::error_info<struct lns_sec_name,std::string> section_name;
      //! Name of the section
      typedef boost::error_info<struct lns_op_name,std::string> option_name;
    }
  }
}

#endif

#ifndef LADA_LNS_EXCEPTIONS_H_
#define LADA_LNS_EXCEPTIONS_H_

#include "LaDaConfig.h"

#include <iostream>

#include <boost/throw_exception.hpp>
#include <boost/exception/error_info.hpp>
#include <boost/exception/info.hpp>


namespace LaDa
{
  namespace crystal
  {
    namespace error
    {
      //! Internal error.
      struct internal: virtual boost::exception, virtual std::exception {};

      //! Root of input errors.
      struct input: virtual boost::exception, virtual std::exception {};
      //! Thrown when a function requires the structure have only one specie per atom.
      struct too_many_species: virtual input {};

//     //! Name of the section
//     typedef boost::error_info<struct lns_sec_name,std::string> section_name;
//     //! Name of the section
//     typedef boost::error_info<struct lns_op_name,std::string> option_name;
    }
  }
}

#endif

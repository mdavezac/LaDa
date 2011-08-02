#ifndef LADA_LNS_EXCEPTIONS_H_
#define LADA_LNS_EXCEPTIONS_H_

#include "LaDaConfig.h"

#include <iostream>

#include <boost/throw_exception.hpp>
#include <boost/exception/error_info.hpp>
#include <boost/exception/info.hpp>


namespace LaDa
{
  namespace math
  {
    namespace error
    {
      //! Base of math exceptions.
      struct math: virtual boost::exception, virtual std::exception {};

      //! Root of input errors.
      struct array_of_different_sizes: virtual math {};
//     //! Name of the section
//     typedef boost::error_info<struct lns_sec_name,std::string> section_name;
//     //! Name of the section
//     typedef boost::error_info<struct lns_op_name,std::string> option_name;
    }
  }
}

#endif

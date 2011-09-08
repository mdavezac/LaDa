#ifndef LADA_ROOTEXCEPTIONS_H
#define LADA_ROOTEXCEPTIONS_H

#include "LaDaConfig.h"

#include <string>
#include <stdexcept>

#include <boost/throw_exception.hpp>
#include <boost/exception/error_info.hpp>
#include <boost/exception/info.hpp>


namespace LaDa
{
  namespace error
  {
    //! Root exception for all lada exceptions.
    struct root: virtual boost::exception, virtual std::exception {};

    //! Root of input errors.
    struct input: virtual root {};

    //! \brief Root of internal error.
    //! \brief These should be programmer errors, rather than something the
    //!        users would see.
    struct internal: virtual root {};
    //! out-of-range error.
    struct out_of_range: virtual internal {};
    //! \brief end-of-iteration error.
    //! \details Should be used to avoid infinit loops. 
    struct stop_iteration: virtual internal {};

    //! Convenience error info type to capture strings.
    typedef boost::error_info<struct string_info,std::string> string;
  }
}

#endif

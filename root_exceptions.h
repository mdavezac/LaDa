#ifndef LADA_ROOTEXCEPTIONS_H
#define LADA_ROOTEXCEPTIONS_H

#include "LaDaConfig.h"

#include <string>
#include <stdexcept>

#include <Python.h>

#include <boost/throw_exception.hpp>
#include <boost/exception/all.hpp>


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
    struct infinite_loop: virtual internal {};

    //! \brief Internal error thrown when encountering a python exception.
    //! \details 
    struct pyerror : virtual internal {};

    //! Convenience error info type to capture strings.
    typedef boost::error_info<struct string_info,std::string> string;
    //! \brief Convenience error infor type to capture python objects.
    //! \details No increffing or decreffing. A python exception should already
    //!          exist.
    typedef boost::error_info<struct string_info,PyObject*> pyobject;
    //! \brief Convenience error infor type to capture python exception.
    //! \details Holds python exception to throw if no PyErr yet.
    typedef boost::error_info<struct string_info,PyObject*> pyexcept;
  }
}

#endif

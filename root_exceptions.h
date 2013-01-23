#ifndef PYLADA_ROOTEXCEPTIONS_H
#define PYLADA_ROOTEXCEPTIONS_H

#include "PyladaConfig.h"

#include <string>
#include <stdexcept>

#include <Python.h>


namespace Pylada
{
  namespace error
  {
    //! Root exception for all pylada exceptions.
    struct root: virtual std::exception {};

    //! Root of input errors.
    struct input: virtual root {};

    //! \brief Root of internal error.
    //! \brief These should be programmer errors, rather than something the
    //!        users would see.
    struct internal: virtual root {};
    //! out-of-range error.
    struct out_of_range: virtual root {};
    //! \brief end-of-iteration error.
    //! \details Should be used to avoid infinit loops. 
    struct infinite_loop: virtual root {};

    //! Convenience error info type to capture strings.
    // typedef boost::error_info<struct string_info,std::string> string;
    //! \brief Convenience error infor type to capture python objects.
    //! \details No increffing or decreffing. A python exception should already
    //!          exist.
    typedef void * pyobject;

  }

  void BOOST_THROW_EXCEPTION( error::root exc);

}

#endif

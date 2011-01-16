#ifndef LADA_VFF_EXCEPTIONS_H_
#define LADA_VFF_EXCEPTIONS_H_

#include "LaDaConfig.h"

#include <iostream>

#include <boost/exception/error_info.hpp>
#include <boost/exception/info.hpp>


namespace LaDa
{
  namespace vff
  {
    //! Internal error.
    struct internal: virtual boost::exception, virtual std::exception {}; 
    //! Throws this when input is somehow incorrect.
    struct input: virtual boost::exception, virtual std::exception {}; 
    //! Throws this when bond input is somehow incorrect.
    struct bond_input: virtual input, virtual std::exception {}; 
    //! Throws this when angle input is somehow incorrect.
    struct angle_input: virtual input, virtual std::exception {}; 
    //! Throws this when wrong number of bonds is encountered.
    struct wrong_structure: virtual input, virtual std::exception {}; 

    //! Qualifies errors.
    typedef boost::error_info<struct vff_error, std::string> error_string; 
  }
}

#endif

#ifndef LADA_VFF_EXCEPTIONS_H_
#define LADA_VFF_EXCEPTIONS_H_

#include "LaDaConfig.h"

#include <iostream>

#include <boost/exception/error_info.hpp>
#include <boost/exception/info.hpp>
#include <boost/tuple/tuple.hpp>

#include <crystal/structure.h>


namespace LaDa
{
  namespace vff
  {
    namespace exceptions
    {
      //! Internal error.
      struct internal: virtual boost::exception, virtual std::exception {}; 
      //! Thrown when input is somehow incorrect.
      struct input: virtual boost::exception, virtual std::exception {}; 
      //! Thrown when bond input is somehow incorrect.
      struct bond_input: virtual input, virtual std::exception {}; 
      //! Thrown when angle input is somehow incorrect.
      struct angle_input: virtual input, virtual std::exception {}; 
      //! Thrown on missing bond parameters.
      struct missing_bond: virtual input, virtual std::exception {};
      //! Thrown on missing angle  parameters.
      struct missing_angle: virtual input, virtual std::exception {};
      //! Thrown this when wrong number of bonds is encountered.
      struct faulty_structure: virtual input, virtual std::exception {}; 
      //! Thrown inside first neighbor tree building.
      struct tree_building: virtual boost::exception, virtual std::exception {};
      //! Thrown inside first neighbor tree building.
      struct site_index: virtual tree_building {};
      //! Thrown when origins don't match in angle_key member function.
      struct mishappened_origin : virtual internal {};
      
      //! Qualifies errors.
      typedef boost::error_info<struct vff_error, std::string> string; 
      //! Qualifies errors.
      typedef boost::error_info<struct vff_error, int> integer; 
      //! Qualifies errors.
      typedef boost::error_info<struct vff_error, boost::tuple<std::string, std::string> > bond; 
      //! Qualifies errors.
      typedef boost::error_info<struct vff_error, boost::tuple<std::string, std::string, std::string> > angle; 
      //! Qualifies errors.
      typedef boost::error_info<struct vff_error, Crystal::TStructure<std::string>::t_Atom > atom; 
    }
  }
}

#endif

//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <sstream>
#include <complex>

#include <boost/exception/diagnostic_information.hpp>

#include <boost/python/class.hpp>
#include <boost/python/errors.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <boost/python/def.hpp>

#define LADA_PYTHON_STD_VECTOR_NOPRINT
#include <python/std_vector.hpp>
#include <crystal/lattice.h>
#include <python/misc.hpp>

#include <crystal/lattice.h>

#include "translation.hpp"
#include "../exceptions.h"
#include "../translation.h"



namespace LaDa
{
  
  namespace Python
  {
    namespace bp = boost::python;
    void expose_translation()
    {
      bp::def
      (
        "create_translations",
        &enumeration::create_translations,
        (
          bp::arg("smith"),
          bp::arg("nsites")
        ),
        "Returns array of translation operators for a given "
        "Smith normal form and number of lattice sites."
      );

      bp::scope scope = bp::class_<enumeration::Translation>
      (
        "Transform", 
        "Rotation + translation.", 
        bp::init<atat::iVector3d const &, atat::iVector3d const&, size_t>()
      ).def( bp::init<enumeration::Translation const&>() )
       .def("__call__", &enumeration::Translation::operator());


      expose_vector<enumeration::Translation>("Array", "Array of Translations");
      bp::register_ptr_to_python< boost::shared_ptr< std::vector<enumeration::Translation> > >();
    }

  }
} // namespace LaDa

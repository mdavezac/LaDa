#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <string>
#include <sstream>

#include <boost/python/class.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/errors.hpp>

#include "../redirect.h"

#include "redirect.hpp"


namespace LaDa
{
  namespace python
  {
    namespace bp = boost::python;
    opt::FortranRedirect *create( opt::FortranRedirect::Unit _unit,
                                  std::string const &_path, 
                                  bool _append )
    {
      opt::FortranRedirect *result = new opt::FortranRedirect(_unit, _path, _append);
      if( not result )
      {
        PyErr_SetString(PyExc_IOError, "Could not redirect fortran logical unit.");
        bp::throw_error_already_set();
        return NULL;
      }
      return result;
    }

    void expose_errors()
    {
      using namespace boost::python;
      bp::scope scope = class_< opt::FortranRedirect, boost::noncopyable >
        ( 
           "Redirect", 
           "Redirects fortran pre-connected units to file.\n\n",
           bp::no_init
        ).def("__init__", bp::make_constructor(&create)  )
         .def("close", &opt::FortranRedirect::close)
         .add_property("is_open", &opt::FortranRedirect::is_open);

      bp::enum_<opt::FortranRedirect::Unit>("fortran", "Fortran pre-connected unit to redirect.\n")
        .value("input", opt::FortranRedirect::input)
        .value("output", opt::FortranRedirect::output)
        .value("error", opt::FortranRedirect::error)
        .export_values();
    }
  }
} // namespace LaDa

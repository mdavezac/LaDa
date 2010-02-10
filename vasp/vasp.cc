#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python/def.hpp>
#include <boost/python/module.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/scope.hpp>

extern "C" { FC_FUNC(vasplib, vasplib)(int *); }

BOOST_PYTHON_MODULE(_vasp)
{
  namespace bp = boost::python;
  bp::scope scope;
  scope.attr("__doc__") = "VASP-as-library extension.";
  bp::docstring_options doc_options(true, false);

  bp::def( "vasp", &FC_FUNC(vasplib, vasplib),
           "Make a call to vasp, given the handle to a communicator." );
}

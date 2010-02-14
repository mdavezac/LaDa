#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <string>
#include <sstream>

#include <boost/python/class.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/default_call_policies.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/errors.hpp>
#include <boost/python/return_arg.hpp>

#include "../redirect.h"

#include "redirect.hpp"


namespace LaDa
{
  namespace python
  {
    namespace bp = boost::python;
    struct Redirect : public boost::noncopyable
    {
      opt::FortranRedirect::Unit unit;
      std::string path;
      bool append;

      Redirect   ( opt::FortranRedirect::Unit _unit, std::string const _path, bool _append)
               : unit(_unit), path(_path), append(_append)
      {
        if( path.empty() ) 
          switch(unit)
          {
            case opt::FortranRedirect::input: path = "stdin"; break;
            case opt::FortranRedirect::output: path = "stdout"; break;
            case opt::FortranRedirect::error: path = "stderr"; break;
          } 
      }
      void __enter__()
      { 
        boost::shared_ptr<opt::FortranRedirect>
          impl(new opt::FortranRedirect(unit, path, append));
        if( not impl )
        {
          PyErr_SetString(PyExc_IOError, "Could not redirect fortran logical unit.");
          bp::throw_error_already_set();
          return;
        }
        impl_.swap(impl);
      }
      void __exit__(bp::object const&, bp::object const&, bp::object const&)
      {
        boost::shared_ptr<opt::FortranRedirect> null_;
        impl_.swap(null_);
      }
      private:
        boost::shared_ptr<opt::FortranRedirect> impl_;
    };

    void expose_redirect()
    {
      using namespace boost::python;
      bp::scope scope = class_< Redirect, boost::noncopyable >
        ( 
           "Redirect", 
           "Redirects fortran pre-connected units to file.\n\n"
           "Should be used with \"with\" context statements:\n\n"
           ">>> with Redirected(Redirected.fortran.ouput, \"stdout\") as redirected:\n"
           ">>>   # make a call to a fortran library in this context.\n\n"
           "@param unit: Fortran pre-connected unit to redirect.\n"
           "@type unit: L{fortran}\n"
           "@param filename: filename to which to redirect.\n"
           "@type filename: str\n"
           "@param append: If true, appends to file, truncates file if false.\n"
           "@type append: boolean\n",
           bp::init< opt::FortranRedirect::Unit const, std::string const&, bool>
              ((bp::arg("unit"), bp::arg("filename"), bp::arg("append")=false))
        ).def("__exit__", &Redirect::__exit__)
         .def("__enter__", &Redirect::__enter__, bp::return_self<>());

      bp::enum_<opt::FortranRedirect::Unit>("fortran", "Fortran pre-connected unit to redirect.\n")
        .value("input", opt::FortranRedirect::input)
        .value("output", opt::FortranRedirect::output)
        .value("error", opt::FortranRedirect::error)
        .export_values();
    }
  }
} // namespace LaDa

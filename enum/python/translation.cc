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
#include <boost/python/return_internal_reference.hpp>


#include <crystal/lattice.h>

#include "translation.hpp"
#include "../exceptions.h"
#include "../translation.h"



namespace LaDa
{
  
  namespace Python
  {
    namespace bp = boost::python;
    struct Translation : enumeration::Translation
    {
      Translation   (atat::iVector3d const &_smith, size_t _nsites)
                  : enumeration::Translation(_smith, _nsites), first_(true) {}
      Translation   (const Translation &_c)
                  : enumeration::Translation(_c), first_(_c.first_) {}

      Translation &iter()  { return *this; }
      Translation const& next()
      {
        if( first_ ) first_ = false; 
        else if( not enumeration::Translation::operator++() )
        {
          first_ = true;
          PyErr_SetString(PyExc_StopIteration, "End of translations.");
          bp::throw_error_already_set();
        }
        return *this;
      }
      bool first_;
    };

    void expose_translation()
    {
      bp::class_<Translation>
      (
        "Translation", 
        "Rotation + translation.", 
        bp::init<atat::iVector3d const &, size_t>()
      ).def( bp::init<Translation const&>() )
       .def("__call__", &Translation::operator())
       .def("__len__", &Translation::size)
       .def("__iter__", &Translation::iter, bp::return_internal_reference<1>())
       .def("next", &Translation::next, bp::return_internal_reference<1>());
    }

  }
} // namespace LaDa

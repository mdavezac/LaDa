//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <sstream>
#include <complex>

#include <boost/python/class.hpp>
#include <boost/python/errors.hpp>
#include <boost/python/return_internal_reference.hpp>

#include <crystal/lattice.h>

#include <crystal/lattice.h>

#include "label_exchange.hpp"
#include "../exceptions.h"
#include "../label_exchange.h"



namespace LaDa
{
  
  namespace Python
  {
    namespace bp = boost::python;

    enumeration::LabelExchange & iter( enumeration::LabelExchange &_self ) { return _self; }
    enumeration::LabelExchange const& next( enumeration::LabelExchange &_self )
    {
      if( not ++_self )
      {
        PyErr_SetString( PyExc_StopIteration, "End of range.\n" );
        bp::throw_error_already_set();
      }
      return _self; 
    }

    void expose_label_exchange()
    {
      bp::class_<enumeration::LabelExchange>
      (
        "LabelExchange", 
        "Rotation + translation.", 
        bp::init<size_t, size_t>()
      ).def( bp::init<enumeration::LabelExchange const&>() )
       .def("__iter__", &iter, bp::return_internal_reference<1>())
       .def("next", &next, bp::return_internal_reference<1>())
       .def("__call__", &enumeration::LabelExchange::operator());
    }

  }
} // namespace LaDa

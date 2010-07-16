#include "LaDaConfig.h"

#include <sstream>
#include <complex>

#include <boost/exception/diagnostic_information.hpp>

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
      if( not (++_self) )
      {
        PyErr_SetString( PyExc_StopIteration, "End of range.\n" );
        bp::throw_error_already_set();
      }
      return _self; 
    }

    inline enumeration::t_uint __call__( enumeration::LabelExchange const &_self, 
                                         enumeration::t_uint _x, 
                                         enumeration::FlavorBase _fl )
    {
      try { return _self(_x, _fl); } 
      catch(enumeration::integer_too_large &_e) 
      {
        std::ostringstream sstr;
        sstr << "Argument integer too large.\n"
             << boost::diagnostic_information(_e);
        PyErr_SetString( PyExc_ValueError, sstr.str().c_str());
        bp::throw_error_already_set();
      }
      catch(enumeration::argument_error &_e) 
      {
        PyErr_SetString( PyExc_ValueError, "Incorrect Argument. ");
        bp::throw_error_already_set();
      }
      catch(std::exception &_e)
      {
        PyErr_SetString( PyExc_RuntimeError, 
                         ("Caught C++ exception: " + std::string(_e.what()) + ".\n").c_str() );
        bp::throw_error_already_set();
      }
      return enumeration::t_uint(-1);
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
       .def("__call__", &__call__);
    }

  }
} // namespace LaDa

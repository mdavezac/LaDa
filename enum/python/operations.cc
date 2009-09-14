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

#include "operations.hpp"
#include "../translation.h"
#include "../label_exchange.h"
#include "../transform.h"



namespace LaDa
{
  
  namespace Python
  {
    namespace bp = boost::python;
    void init( enumeration::Transform &_self, bp::tuple const &_tuple )
    {
      if( bp::len(_tuple) != 2 )
      {
        PyErr_SetString(PyExc_TypeError, "Argument is not a 2-tuple.\n");
        bp::throw_error_already_set();
        return;
      }
      try
      {
        Crystal::t_SmithTransform t;
        boost::tuples::get<0>(t) =  bp::extract<atat::rMatrix3d>( _tuple[0] );
        boost::tuples::get<1>(t) =  bp::extract<atat::iVector3d>( _tuple[1] );
        _self.init(t); 
      }
      catch(boost::exception &_e)
      {
        std::ostringstream sstr;
        sstr << "Internal error found while initializing enumeration.Transform object.\n"
             << boost::diagnostic_information(_e);
        PyErr_SetString(PyExc_RuntimeError, sstr.str().c_str());
        bp::throw_error_already_set();
        return;
      }
      catch(...)
      {
        PyErr_SetString(PyExc_RuntimeError, "Could not extract tuple.\n");
        bp::throw_error_already_set();
        return;
      }
    }

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

    void expose_operations()
    {
      bp::class_<enumeration::Transform>
      (
        "Transform", 
        "Rotation + translation.", 
        bp::init<Crystal::SymmetryOperator const &, Crystal::Lattice const &>()
      ).def( bp::init<enumeration::Transform const&>() )
       .def("init", &init)
       .def("__call__", &enumeration::Transform::operator());

      bp::class_<enumeration::LabelExchange>
      (
        "Transform", 
        "Rotation + translation.", 
        bp::init<size_t, size_t>()
      ).def( bp::init<enumeration::LabelExchange const&>() )
       .def("__iter__", &iter, bp::return_internal_reference<1>())
       .def("next", &next, bp::return_internal_reference<1>())
       .def("__call__", &enumeration::LabelExchange::operator());

      bp::def
      (
        "create_translation",
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

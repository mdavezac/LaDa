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
#include "../exceptions.h"
#include "../translation.h"
#include "../label_exchange.h"
#include "../transform.h"



namespace LaDa
{
  
  namespace Python
  {
    namespace bp = boost::python;
    boost::shared_ptr< std::vector<enumeration::Transform> > 
      create_transforms( Crystal::Lattice const & _lat)
      {
        try { return enumeration::create_transforms(_lat); }
        catch(enumeration::symmetry_not_of_lattice &_e)
        {
          std::ostringstream sstr;
          sstr << "Rotation+Translation does not leave the lattice invariant.\n"
               << boost::diagnostic_information(_e);
          PyErr_SetString(PyExc_RuntimeError, sstr.str().c_str());
          bp::throw_error_already_set();
          return boost::shared_ptr< std::vector<enumeration::Transform> >();
        }
        catch(boost::exception &_e)
        {
          std::ostringstream sstr;
          sstr << "Internal error found while creating array of Transforms.\n"
               << boost::diagnostic_information(_e);
          PyErr_SetString(PyExc_RuntimeError, sstr.str().c_str());
          bp::throw_error_already_set();
          return boost::shared_ptr< std::vector<enumeration::Transform> >();
        }
        catch(...)
        {
          PyErr_SetString(PyExc_RuntimeError, "Uknown error.\n");
          bp::throw_error_already_set();
          return boost::shared_ptr< std::vector<enumeration::Transform> >();
        }
      }

    void init( enumeration::Transform &_self,
               atat::rMatrix3d const& _left,
               atat::iVector3d const& _smith)
    {
      try { _self.init(_left, _smith); }
      catch(enumeration::symmetry_not_of_supercell &_e)
      {
        std::ostringstream sstr;
        sstr << "Rotation+Translation does not leave the supercell invariant.\n"
             << boost::diagnostic_information(_e);
        PyErr_SetString(PyExc_RuntimeError, sstr.str().c_str());
        bp::throw_error_already_set();
        return;
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
        PyErr_SetString(PyExc_RuntimeError, "Unkown error.");
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

    void expose_transform()
    {
      bp::def
      (
        "create_transforms",
        &create_transforms,
        bp::arg("lattice"),
        "Returns an array of rotation+translation operators for a given Lattice."
      );

      bp::class_<enumeration::Transform>
      (
        "Transform", 
        "Rotation + translation. Must be initialized to appropriate Smith normal form.", 
        bp::init<Crystal::SymmetryOperator const &, Crystal::Lattice const &>()
      ).def( bp::init<enumeration::Transform const&>() )
       .def("init", &init, (bp::arg("left"), bp::arg("smith")),
            "Initializes the transform for a specific supercell.\n"
            " _ left is the left transform matrix from the Hermite to the Smith normal form.\n"
            " _ smith is the diagonal of the Smith normal form as an atat.iVector3d.\n")
       .def("__call__", &enumeration::Transform::operator())
       .def_readwrite("op", &Crystal::SymmetryOperator::op)
       .def_readwrite("trans", &Crystal::SymmetryOperator::trans)
       .def("invariant", &Crystal::SymmetryOperator::invariant, 
            (bp::arg("matrix"), bp::arg("tolerance")=types::tolerance),
            "Returns true if the matrix is invariant through this rotation.")
       .def("__call__", &Crystal::SymmetryOperator::SymmetryOperator::operator())
       .def("__str__", &tostream<Crystal::SymmetryOperator>);

      expose_vector<enumeration::Transform>("Array", "Array of Translations");
      bp::register_ptr_to_python< boost::shared_ptr< std::vector<enumeration::Transform> > >();
    }
  }
} // namespace LaDa

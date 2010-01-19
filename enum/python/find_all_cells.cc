//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <sstream>
#include <complex>

#include <boost/python/class.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/errors.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

#include <python/std_vector.hpp>
#include <crystal/lattice.h>
#include <python/misc.hpp>

#include "find_all_cells.hpp"
#include "../find_all_cells.h"



namespace LaDa
{
  
  namespace Python
  {
    boost::shared_ptr< std::vector<enumeration::SmithGroup> >
      create_smith_groups1( Crystal::Lattice const &_lattice, size_t _nmax)
        { return enumeration::create_smith_groups( _lattice, _nmax); }
    boost::shared_ptr< std::vector<enumeration::SmithGroup> >
      create_smith_groups2( Crystal::Lattice const &_lattice,
                            boost::shared_ptr< std::vector<math::rMatrix3d> > const & _s )
        { return enumeration::create_smith_groups( _lattice, _s); }

    enumeration::SmithGroup::Supercell* create( boost::python::tuple const &_tuple)
    {
      namespace bp = boost::python;
      if( bp::len(_tuple) != 2 )
      {
        PyErr_SetString(PyExc_TypeError, "Argument is not a 2-tuple.\n");
        bp::throw_error_already_set();
        return new enumeration::SmithGroup::Supercell();
      }
      try
      {
        math::rMatrix3d const transform = bp::extract<math::rMatrix3d>( _tuple[0] );
        math::rMatrix3d const hermite = bp::extract<math::rMatrix3d>( _tuple[1] );
        return new enumeration::SmithGroup::Supercell(transform, hermite);
      }
      catch(...)
      {
        PyErr_SetString(PyExc_RuntimeError, "Could not extract tuple.\n");
        bp::throw_error_already_set();
        return new enumeration::SmithGroup::Supercell();
      }
    }
    void expose_find_all_cells()
    {
      namespace bp = boost::python;
      bp::def
      (
        "find_all_cells", 
        &LaDa::enumeration::find_all_cells,
        "Finds all cells of a certain size for a given lattice."
      );

      bp::register_ptr_to_python< boost::shared_ptr< std::vector<math::rMatrix3d> > >();

      bp::def
      (
        "create_smith_groups",
        &create_smith_groups1, 
        "Finds all inequivalent matrices of size N, grouped per Smith normal transform.\n",
        (
          bp::arg("lattice"),
          bp::arg("N")
        )
      );
      bp::def
      (
        "create_smith_groups",
        &create_smith_groups2, 
        "Groups hermite matrices per Smith normal transform.",
        (
          bp::arg("lattice"),
          bp::arg("hermites")
        )
      );

      typedef enumeration::SmithGroup t_SG;
      bp::scope scope = bp::class_<enumeration::SmithGroup>
      (
        "SmithGroup", 
        "A group of supercells with equivalent translational symmetries",
        bp::init<math::iVector3d const&>()
      ).def(bp::init<enumeration::SmithGroup const &>())
       .add_property
       (
         "smith",
         make_getter(&t_SG::smith, bp::return_value_policy<bp::return_by_value>()),
         make_setter(&t_SG::smith, bp::return_value_policy<bp::return_by_value>()),
         "Smith translational group.\n\nNumpy integer 3x1 array."
       ) 
       .def_readwrite("supercells", &enumeration::SmithGroup::supercells)
       .def("__str__", &tostream<enumeration::SmithGroup>);

      bp::class_<enumeration::SmithGroup::Supercell>
      (
        "Supercell", 
        "A supercell defined by a Hermite matrix and a transform to translation group.\n",
        bp::init<math::rMatrix3d const&, math::rMatrix3d const&>()
      ).def(bp::init<enumeration::SmithGroup::Supercell const&>())
       .def( "__init__", bp::make_constructor( &create ) )
       .add_property
       (
         "transform",
         make_getter(&t_SG::Supercell::transform, bp::return_value_policy<bp::return_by_value>()),
         make_setter(&t_SG::Supercell::transform, bp::return_value_policy<bp::return_by_value>()),
         "Transform matrix between Smith representation and hermite representation.\n\n"
         "Numpy float64 3x3 array."
       ) 
       .add_property
       (
         "hermite",
         make_getter(&t_SG::Supercell::hermite, bp::return_value_policy<bp::return_by_value>()),
         make_setter(&t_SG::Supercell::hermite, bp::return_value_policy<bp::return_by_value>()),
         "Hermite matrix.\n\nNumpy float64 3x3 array."
       ) 
       .def("__str__", &tostream<enumeration::SmithGroup::Supercell>);

      
      expose_vector<enumeration::SmithGroup::Supercell>
         ("SupercellsArray", "An array of enumeration.SmithGroup");
      expose_vector<enumeration::SmithGroup>
         ("Array", "An array of enumeration.SmithGroup");
      bp::register_ptr_to_python< boost::shared_ptr< std::vector<enumeration::SmithGroup> > >();
     

    }

  }
} // namespace LaDa

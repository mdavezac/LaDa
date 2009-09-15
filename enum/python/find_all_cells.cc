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
                            boost::shared_ptr< std::vector<atat::rMatrix3d> > const & _s )
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
        atat::rMatrix3d const transform = bp::extract<atat::rMatrix3d>( _tuple[0] );
        atat::rMatrix3d const hermite = bp::extract<atat::rMatrix3d>( _tuple[1] );
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
      expose_vector<atat::rMatrix3d>
         ("rMatrix3dArray", "An array of atat.rMatrix3dArray");

      bp::def
      (
        "find_all_cells", 
        &LaDa::enumeration::find_all_cells,
        "Finds all cells of a certain size for a given lattice."
      );

      bp::register_ptr_to_python< boost::shared_ptr< std::vector<atat::rMatrix3d> > >();

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

      bp::scope scope = bp::class_<enumeration::SmithGroup>
      (
        "SmithGroup", 
        "A group of supercells with equivalent translational symmetries",
        bp::init<atat::iVector3d const&>()
      ).def(bp::init<enumeration::SmithGroup const &>())
       .def_readwrite("smith", &enumeration::SmithGroup::smith)
       .def_readwrite("supercells", &enumeration::SmithGroup::supercells)
       .def("__str__", &tostream<enumeration::SmithGroup>);

      bp::class_<enumeration::SmithGroup::Supercell>
      (
        "Supercell", 
        "A supercell defined by a Hermite matrix and a transform to translation group.\n",
        bp::init<atat::rMatrix3d const&, atat::rMatrix3d const&>()
      ).def(bp::init<enumeration::SmithGroup::Supercell const&>())
       .def( "__init__", bp::make_constructor( &create ) )
       .def_readwrite("transform", &enumeration::SmithGroup::Supercell::transform)
       .def_readwrite("hermite", &enumeration::SmithGroup::Supercell::hermite)
       .def("__str__", &tostream<enumeration::SmithGroup::Supercell>);

      
      expose_vector<enumeration::SmithGroup::Supercell>
         ("SupercellsArray", "An array of enumeration.SmithGroup");
      expose_vector<enumeration::SmithGroup>
         ("Array", "An array of enumeration.SmithGroup");
      bp::register_ptr_to_python< boost::shared_ptr< std::vector<enumeration::SmithGroup> > >();
     

    }

  }
} // namespace LaDa

#include "LaDaConfig.h"

#include <sstream>
#include <complex>

#include <boost/python.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/def.hpp>
#include <boost/python/errors.hpp>

#include "../smith.h"


namespace LaDa
{
  namespace Python
  {

    
    boost::python::tuple get_smith_transform( math::rMatrix3d const &_lat_cell,
                                              math::rMatrix3d const &_str_cell )
    {
      namespace bt = boost::tuples;
      namespace bp = boost::python;
      try
      {
        Crystal::t_SmithTransform result = Crystal::get_smith_transform( _lat_cell, _str_cell );
        return bp::make_tuple( bt::get<0>(result), bt::get<1>(result) );
      }
      catch(...)
      {
        PyErr_SetString
        (
          PyExc_IOError, 
          "Error while computing transform to smith normal form.\n" 
        );
        bp::throw_error_already_set();
        return bp::make_tuple( -1 );
      }
    }
    template< class T_TYPE  >
    boost::python::tuple get_smith_transform_str( Crystal::TStructure<T_TYPE> const &_struct )
    {
      namespace bt = boost::tuples;
      namespace bp = boost::python;
      try
      {
        Crystal::t_SmithTransform result = Crystal::get_smith_transform( _struct );
        return bp::make_tuple( bt::get<0>(result), bt::get<1>(result) );
      }
      catch(...)
      {
        PyErr_SetString
        (
          PyExc_IOError, 
          "Error while computing transform to smith normal form.\n" 
        );
        bp::throw_error_already_set();
        return bp::make_tuple( -1 );
      }
    }

    
    //! Computes smith indices of position \a _pos.
    boost::python::tuple get_smith_index( boost::python::tuple const & _transform,
                                          math::rVector3d  const &_pos )
    {
      namespace bt = boost::tuples;
      namespace bp = boost::python;
      if( bp::len( _transform ) != 2 )
      {
        PyErr_SetString
        (
          PyExc_RuntimeError, 
          "Incorrect tuple argument when computing smith normal form indices.\n" 
        );
        bp::throw_error_already_set();
        return bp::make_tuple(-1,-1,-1);
      }
      try
      {
        Crystal::t_SmithTransform transform;
        try{  bt::get<0>( transform ) = bp::extract< math::rMatrix3d >( _transform[0] ); }
        catch(...)
        {
          PyErr_SetString
          (
            PyExc_IOError, 
            "First tuple element is not an rMatrix3d object.\n"
          );
          bp::throw_error_already_set();
          return bp::make_tuple(-1,-1,-1);
        }
        try{  bt::get<1>( transform ) = bp::extract< math::iVector3d >( _transform[1] ); }
        catch(...)
        {
          PyErr_SetString
          (
            PyExc_IOError, 
            "First tuple element is not an rMatrix3d object.\n"
          );
          bp::throw_error_already_set();
          return bp::make_tuple(-1,-1,-1);
        }
        math::iVector3d const vec( Crystal::get_smith_index( transform, _pos ) );
        return bp::make_tuple( vec(0), vec(1), vec(2) );
      }
      catch(...)
      {
        PyErr_SetString
        (
          PyExc_IOError, 
          "Error while computing smith normal form indices.\n" 
        );
        bp::throw_error_already_set();
        return bp::make_tuple(-1,-1,-1);
      }
    }

    //! Computes smith indices of position \a _pos.
    int get_linear_smith_index( boost::python::tuple const & _transform,
                                math::rVector3d  const &_pos, int _site_index )
    {
      namespace bt = boost::tuples;
      namespace bp = boost::python;
      if( bp::len( _transform ) != 2 )
      {
        PyErr_SetString
        (
          PyExc_TypeError, 
          "Incorrect tuple argument when computing smith normal form indices.\n" 
        );
        bp::throw_error_already_set();
        return -1;
      }
      if( _site_index < 0 )
      {
        PyErr_SetString
        (
          PyExc_ValueError, 
          "Negative site index not accepted in linear_smith_index."
        );
        bp::throw_error_already_set();
        return -1;
      }
      Crystal::t_SmithTransform transform;
      try{  bt::get<0>( transform ) = bp::extract< math::rMatrix3d >( _transform[0] ); }
      catch(...)
      {
        PyErr_SetString
        (
          PyExc_TypeError, 
          "First tuple element is not an rMatrix3d object.\n"
        );
        bp::throw_error_already_set();
        return -1;
      }
      try{  bt::get<1>( transform ) = bp::extract< math::iVector3d >( _transform[1] ); }
      catch(...)
      {
        PyErr_SetString
        (
          PyExc_TypeError, 
          "Second tuple element is not an iVector3d object.\n"
        );
        bp::throw_error_already_set();
        return -1;
      }
      try{ return Crystal::get_linear_smith_index( transform, _site_index, _pos ); }
      catch(...)
      {
        PyErr_SetString
        (
          PyExc_RuntimeError, 
          "Error while computing smith normal form indices.\n" 
        );
        bp::throw_error_already_set();
        return -1;
      }
    }

    void expose_smith()
    {
      namespace bp = boost::python;
      bp::def( "smith_normal_transform", &get_smith_transform_str<std::string>,
               bp::arg("structure") );
      bp::def( "smith_normal_transform", &get_smith_transform_str<types::t_real>,
               bp::arg("structure") );
      bp::def
      ( 
        "smith_normal_transform", &get_smith_transform,
        (bp::arg("unitcell"), bp::arg("supercell")), 
        "Computes Smith Normal Transform parameters.\n\n" 
	      "  structure : `Structure`\n    A structure for which *structure.lattice* is set. "
          "The supercell is then *structure.cell* and the unit cell "
          "*structure.lattice.cell*. \n"
        "  unitcell\n    Primitive cell of the lattice from which "
           "supercells are built.\n"
        "  supercell\n    Multiple of the unit cell.\n\n"
        ":return: tuple allowing a stransformation to the smith normal form.\n\n" 
	      "This function accepts two types of call: using a correctly defined "
          "*structure* alone, and by specifying the *unitcell* and the "
	        "*supercell* explicitely. In any case, the cell of structure must be "
          "exactly commensurate with the lattice, eg no relaxation.\n"
        ":see: `G. Hart and R. Forcade, Phys. Rev. B. **80**, 014120 (2009)"
          "<https:\\dx.doi.org/10.1103/PhysRevB.80.014120>`__\n"
      );
      bp::def
      ( 
        "smith_indices", &get_smith_index,
        ( bp::arg("transform"), bp::arg("position") ),
        "Returns the indices of the position in the smith normal form.\n\n" 
        ":Parameters:"
        "  transform\n    transformation tuple yielded by `smith_normal_transform`\n"
        "  position\n    cartesian coordinates on the lattice.\n"
      );
      bp::def
      ( 
        "linear_smith_index", &get_linear_smith_index,
        ( bp::arg("transform"), bp::arg("position"), bp::arg("site_index") = 0 ),
        "Returns the indices of the position in the smith normal form." 
      );
    }

  } // namespace Python
} // namespace LaDa

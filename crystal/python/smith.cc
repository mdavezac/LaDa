//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

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

    
    boost::python::tuple get_smith_transform( atat::rMatrix3d const &_lat_cell,
                                              atat::rMatrix3d const &_str_cell )
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

    
    //! Computes smith indices of position \a _pos.
    boost::python::tuple get_smith_index( boost::python::tuple const & _transform,
                                          atat::rVector3d  const &_pos )
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
        try{  bt::get<0>( transform ) = bp::extract< atat::rMatrix3d >( _transform[0] ); }
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
        try{  bt::get<1>( transform ) = bp::extract< atat::iVector3d >( _transform[1] ); }
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
        atat::iVector3d const vec( Crystal::get_smith_index( transform, _pos ) );
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

    void expose_smith()
    {
      namespace bp = boost::python;
      bp::def
      ( 
        "smith_normal_transform", &get_smith_transform,
        ( bp::arg("lattice_cell"), bp::arg("structure_cell") ),
        "Returns a tuple allowing a stransformation to the smith normal form." 
      );
      bp::def
      ( 
        "smith_indices", &get_smith_index,
        ( bp::arg("transform"), bp::arg("position") ),
        "Returns the indices of the position in the smith normal form." 
      );
    }

  } // namespace Python
} // namespace LaDa

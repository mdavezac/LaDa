#ifndef LADA_MATH_PYTHON_HPP
#define LADA_MATH_PYTHON_HPP

#include "LaDaConfig.h"

#include <boost/python/object.hpp>
#include <boost/python/extract.hpp>
#include "../eigen.h"

namespace LaDa
{
  namespace python
  {
     //! Checks whether object is convertible to atomic position.
     inline bool is_position(boost::python::object const &_index)
     {
       namespace bp = boost::python;
       if(bp::len(_index) != 3) return false;
       if(not bp::extract<math::rMatrix3d::Scalar>(_index[0]).check()) return false;
       if(not bp::extract<math::rMatrix3d::Scalar>(_index[1]).check()) return false;
       if(not bp::extract<math::rMatrix3d::Scalar>(_index[2]).check()) return false;
       return false;
     }
     //! Converts object to atom. No type checking.
     inline void extract_position(boost::python::object const &_pos, math::rVector3d &_out)
     {
       namespace bp = boost::python;
       _out(0) = bp::extract<math::rMatrix3d::Scalar>(_pos[0]);
       _out(1) = bp::extract<math::rMatrix3d::Scalar>(_pos[0]);
       _out(2) = bp::extract<math::rMatrix3d::Scalar>(_pos[0]);
     }
   }
} 

#endif

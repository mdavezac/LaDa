#include "LaDaConfig.h"

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL math_ARRAY_API
#include <numpy/arrayobject.h>

#include <iostream>

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/object.hpp>
#include <boost/python/errors.hpp>
#include <boost/python/extract.hpp>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Geometry>

#include <python/numpy_types.h>
#include <python/wrap_numpy.h>

#include "avogadro.hpp"
#include "matrix.hpp"
#include "../eigen.h"
#include "../misc.h"

void check_matrix_from_python(boost::python::object const &_object)
{
  try
  {
    LaDa::math::rMatrix3d mat = boost::python::extract<LaDa::math::rMatrix3d>(_object);
    std::cout << mat << "\n";
  }
  catch(std::exception &e)
  {
    std::string error("Caught an error.\n");
    error += e.what();
    PyErr_SetString(PyExc_RuntimeError, error.c_str());
    boost::python::throw_error_already_set();
  }
}

void check_from_python_mat(LaDa::math::rMatrix3d const &_o) { std::cout << _o << "\n"; }
void check_from_python_vec(LaDa::math::rVector3d const &_o) { std::cout << _o << "\n"; }

LaDa::math::rMatrix3d check_topy_matrix()
{
  return LaDa::math::rMatrix3d::Identity();
}
LaDa::math::rVector3d check_topy_vector()
{
  return LaDa::math::rVector3d::Identity();
}
bool is_integer1(LaDa::math::rVector3d const &_in) { return LaDa::math::is_integer(_in); }
bool is_integer2(LaDa::math::rMatrix3d const &_in) { return LaDa::math::is_integer(_in); }

boost::python::object Rotation1(LaDa::types::t_real _angle, boost::python::object const &_a)
{
# ifndef LADA_WITH_EIGEN3 
    // \typedef type of the affine transformations.
    typedef Eigen::Transform<LaDa::types::t_real, 3> Affine;
# else
    // \typedef type of the affine transformations.
    typedef Eigen::Transform<LaDa::types::t_real, 3, Eigen::Isometry> Affine;
# endif
  // \typedef type of the angle axis object to initialize roations.
  typedef Eigen::AngleAxis<LaDa::types::t_real> AngleAxis;
  npy_intp dims[2] = {4, 3};
  int const type = LaDa::math::numpy::type<LaDa::types::t_real>::value;
  PyArrayObject *result = (PyArrayObject*)PyArray_SimpleNew(2, dims, type);
  if(not result) return boost::python::object();
  
  LaDa::math::rVector3d axis;
  LaDa::python::convert_to_vector(_a.ptr(), axis);
  Affine a;
  a = AngleAxis(_angle, axis);
  for(size_t i(0); i < 3; ++i)
    for(size_t j(0); j < 3; ++j)
      *((LaDa::types::t_real*)(result->data + i*result->strides[0] + j*result->strides[1])) = a(i, j);
  for(size_t j(0); j < 3; ++j)
    *((LaDa::types::t_real*)(result->data + 4*result->strides[0] + j*result->strides[1])) = 0;
  return boost::python::object(boost::python::handle<>((PyObject*)result));
}
boost::python::object Translation(boost::python::object const &_a)
{
# ifndef LADA_WITH_EIGEN3 
    // \typedef type of the affine transformations.
    typedef Eigen::Transform<LaDa::types::t_real, 3> Affine;
# else
    // \typedef type of the affine transformations.
    typedef Eigen::Transform<LaDa::types::t_real, 3, Eigen::Isometry> Affine;
# endif
  // \typedef type of the angle axis object to initialize roations.
  npy_intp dims[2] = {4, 3};
  int const type = LaDa::math::numpy::type<LaDa::types::t_real>::value;
  PyArrayObject *result = (PyArrayObject*)PyArray_SimpleNew(2, dims, type);
  if(not result) return boost::python::object();
  
  LaDa::math::rVector3d trans;
  LaDa::python::convert_to_vector(_a.ptr(), trans);
  for(size_t i(0); i < 3; ++i)
    for(size_t j(0); j < 3; ++j)
      *((LaDa::types::t_real*)(result->data + i*result->strides[0] + j*result->strides[1])) = 0;
  for(size_t j(0); j < 3; ++j)
    *((LaDa::types::t_real*)(result->data + 4*result->strides[0] + j*result->strides[1])) = trans(j);
  return boost::python::object(boost::python::handle<>((PyObject*)result));
}


BOOST_PYTHON_MODULE(math)
{
  import_array();
  namespace bp = boost::python;
  bp::scope scope;
  scope.attr("__doc__") = "This namespace defines mathematical types\n\nIt "
   "imports automatic conversions between numpy arrays and 3d-vectors and "
   "matrices.";
  bp::docstring_options doc_options(true, false);

  LaDa::python::expose_eigen_vectors();
  LaDa::python::expose_eigen_matrices();
  bp::def("is_integer", &is_integer1);
  bp::def("is_integer", &is_integer2, "Returns True if the input vector or matrix is integer.");
  bp::def("Rotation", &Rotation1, "Returns rotation according to angle and axis as a 4x3 symmetry operation.");
  bp::def("Translation", &Translation, "Returns translation as last row of a 4x3 symmetry operation.");

  bp::def("check_extract", &check_matrix_from_python, "Check matrix extraction.");
  bp::def("check_frompy_to_cpp_mat", &check_from_python_mat, "Check automatic matrix extraction.");
  bp::def("check_frompy_to_cpp_vec", &check_from_python_vec, "Check automatic vector extraction.");
  bp::def("check_topy_from_cpp_mat", &check_topy_matrix, "Check automatic vector conversion.");
  bp::def("check_topy_from_cpp_vec", &check_topy_vector, "Check automatic vector conversion.");
}

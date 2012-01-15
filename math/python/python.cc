#include "LaDaConfig.h"

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

BOOST_PYTHON_MODULE(math)
{
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

  bp::def("check_extract", &check_matrix_from_python, "Check matrix extraction.");
  bp::def("check_frompy_to_cpp_mat", &check_from_python_mat, "Check automatic matrix extraction.");
  bp::def("check_frompy_to_cpp_vec", &check_from_python_vec, "Check automatic vector extraction.");
  bp::def("check_topy_from_cpp_mat", &check_topy_matrix, "Check automatic vector conversion.");
  bp::def("check_topy_from_cpp_vec", &check_topy_vector, "Check automatic vector conversion.");
}

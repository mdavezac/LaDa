// Last update: timvdm 19 June 2009
// Adapted from Avogadro open source molecular editor from the kde project. (libavogadro/src/python/eigen.cpp)
#include <boost/python/detail/wrap_python.hpp>
#include <numpy/arrayobject.h> 
#include <boost/python.hpp>
#include <boost/python/tuple.hpp>

#include <Eigen/Geometry>
#include "../eigen.h"
#include "avogadro.hpp"

#include <iostream>

namespace LaDa { namespace Python {
using namespace boost::python;

template <typename Scalar> struct ScalarTraits;
template <> struct ScalarTraits<int>
{
  enum { isInt = 1, isFloat = 0, isDouble = 0 };
};
template <> struct ScalarTraits<float>
{
  enum { isInt = 0, isFloat = 1, isDouble = 0 };
};
template <> struct ScalarTraits<double>
{
  enum { isInt = 0, isFloat = 0, isDouble = 1 };
};


  /***********************************************************************
   *
   * Vector3x = Vector3d, Vector3f, Vector3i
   *
   ***********************************************************************/
  
  template <class Vector3x>
  struct Vector3x_to_python_array
  {
    typedef typename Vector3x::Scalar Scalar;
    
    struct innerclass
    {
      //
      //  Eigen::Vector3x --> python array
      //
      static PyObject* convert(Vector3x const &vec)
      {
        npy_intp dims[1] = { 3 };
        PyObject *result;
        if (ScalarTraits<Scalar>::isInt)
          result = PyArray_SimpleNew(1, dims, NPY_INT);
        else if (ScalarTraits<Scalar>::isFloat)
          result = PyArray_SimpleNew(1, dims, NPY_FLOAT);
        else
          result = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
        
        // copy the data
        Scalar *data = (Scalar*) reinterpret_cast<PyArrayObject*>(result)->data;
        data[0] = vec.x();
        data[1] = vec.y();
        data[2] = vec.z();
        
        return incref(result);
      }
 
      //
      //  Eigen::Vector3x * --> python array
      //    
      static PyObject* convert(Vector3x *vec)
      {
        if (!vec)
          throw_error_already_set();
 
        npy_intp dims[1] = { 3 };
        PyObject *result;
        if (ScalarTraits<Scalar>::isInt)
          result = PyArray_SimpleNew(1, dims, NPY_INT);
        else if (ScalarTraits<Scalar>::isFloat)
          result = PyArray_SimpleNew(1, dims, NPY_FLOAT);
        else
          result = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
        
        // copy the data
        Scalar *data = (Scalar*) reinterpret_cast<PyArrayObject*>(result)->data;
        data[0] = vec->x();
        data[1] = vec->y();
        data[2] = vec->z();
        
        return incref(result);
      }

      //
      //  const Eigen::Vector3x * --> python array
      //
      static PyObject* convert(const Vector3x *vec)
      {
        if (!vec)
          throw_error_already_set();

        npy_intp dims[1] = { 3 };
        PyObject *result;
        if (ScalarTraits<Scalar>::isInt)
          result = PyArray_SimpleNew(1, dims, NPY_INT);
        else if (ScalarTraits<Scalar>::isFloat)
          result = PyArray_SimpleNew(1, dims, NPY_FLOAT);
        else
          result = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
        
        // copy the data
        Scalar *data = (Scalar*) reinterpret_cast<PyArrayObject*>(result)->data;
        data[0] = vec->x();
        data[1] = vec->y();
        data[2] = vec->z();
        
        return incref(result);
      }
    };
   
    Vector3x_to_python_array()
    {
      to_python_converter<Vector3x, innerclass>();
      to_python_converter<Vector3x*, innerclass>();
      to_python_converter<const Vector3x*, innerclass>();
    }

  };

  template <class Vector3x>
  struct Vector3x_from_python_array
  {
    typedef typename Vector3x::Scalar Scalar;
    
    Vector3x_from_python_array()
    {
      // Insert an rvalue from_python converter at the tail of the
      // chain. Used for implicit conversions
      //
      //  python array --> Vector3x
      //
      // used for:
      //
      //  void function(Vector3x vec)
      //  void function(Vector3x & vec)
      //  void function(const Vector3x & vec)
      //
      converter::registry::push_back( &convertible, &construct, type_id<Vector3x>() );

      
      // Insert an lvalue from_python converter
      //
      //  python array --> Vector3x*
      //
      // used for:
      //  
      //  void function(const Vector3x * vec)
      converter::registry::insert( &convert, type_id<Vector3x>() );
    }

    static void* convert(PyObject *obj_ptr)
    {
      std::cout << "am here\n";
      if (!PyArray_Check(obj_ptr)) throw_error_already_set();

      // only accept int, long, float and double
      switch (PyArray_ObjectType(obj_ptr, 0))
      {
        case NPY_INT:
        case NPY_LONG:
        case NPY_FLOAT:
        case NPY_DOUBLE:  break;
        default: PyErr_SetString(PyExc_ValueError, "Numpy array has wrong type.\n");
                 throw_error_already_set();
                 return 0;
      }

      // do some type checking
      if ((PyArray_ObjectType(obj_ptr, 0) == NPY_FLOAT) || (PyArray_ObjectType(obj_ptr, 0) == NPY_DOUBLE))
        if (ScalarTraits<Scalar>::isInt)
          return 0;

      if ((PyArray_ObjectType(obj_ptr, 0) == NPY_INT) || (PyArray_ObjectType(obj_ptr, 0) == NPY_LONG))
        if (ScalarTraits<Scalar>::isFloat || ScalarTraits<Scalar>::isDouble)
          return 0;
      
      PyArrayObject *array = reinterpret_cast<PyArrayObject*>(obj_ptr);

      // check the dimensions
      if (array->nd != 1)
      {
        PyErr_SetString(PyExc_ValueError, "Numpy array is not 1d.\n");
        throw_error_already_set(); // the array has at least two dimensions (matrix)
      }
  
      if (array->dimensions[0] != 3)
      {
        PyErr_SetString(PyExc_ValueError, "Numpy array is not of length 3.\n");
        throw_error_already_set(); // the 1D array does not have exactly 3 elements
      }

      int u = PyArray_ObjectType(obj_ptr, 0);
      if( u == NPY_INT )
      {
        int *values = reinterpret_cast<int*>(array->data);
        return new Vector3x(values[0], values[1], values[2]);
      }
      if( u == NPY_LONG )
      {
        long *values = reinterpret_cast<long*>(array->data);
        return new Vector3x(values[0], values[1], values[2]);      
      }
      if( u == NPY_FLOAT )
      {
        float *values = reinterpret_cast<float*>(array->data);
        return new Vector3x(values[0], values[1], values[2]);
      }
      if( u == NPY_DOUBLE )
      {
        double *values = reinterpret_cast<double*>(array->data);
        return new Vector3x(values[0], values[1], values[2]);
      }
      return NULL;
    }
 
    static void* convertible(PyObject *obj_ptr)
    {
      if (!PyArray_Check(obj_ptr)) return 0;

      // only accept int, long, float and double
      switch (PyArray_ObjectType(obj_ptr, 0)) {
        case NPY_INT:
        case NPY_LONG:
        case NPY_FLOAT:
        case NPY_DOUBLE:
          break;
        default:
          return 0;
      }
      
      // do some type checking
      if ((PyArray_ObjectType(obj_ptr, 0) == NPY_FLOAT) || (PyArray_ObjectType(obj_ptr, 0) == NPY_DOUBLE)) 
        if (ScalarTraits<Scalar>::isInt) return 0;

      if ((PyArray_ObjectType(obj_ptr, 0) == NPY_INT) || (PyArray_ObjectType(obj_ptr, 0) == NPY_LONG))
        if (ScalarTraits<Scalar>::isFloat || ScalarTraits<Scalar>::isDouble) return 0;
      
      PyArrayObject *array = reinterpret_cast<PyArrayObject*>(obj_ptr);

      // check the dimensions
      if (array->nd != 1 and array->dimensions[0] != 3) return 0; 
 
      return obj_ptr;
    }

    static void construct(PyObject *obj_ptr, converter::rvalue_from_python_stage1_data *data)
    {
      PyArrayObject *array = reinterpret_cast<PyArrayObject*>(obj_ptr);
      void *storage = ((converter::rvalue_from_python_storage<Vector3x>*)data)->storage.bytes;

      switch (PyArray_ObjectType(obj_ptr, 0)) {
        case NPY_INT:
          {
            int *values = reinterpret_cast<int*>(array->data);
            new (storage) Vector3x(values[0], values[1], values[2]);
          }
          break;
        case NPY_LONG:
          {
            long *values = reinterpret_cast<long*>(array->data);
            new (storage) Vector3x(values[0], values[1], values[2]);
          }
          break;
        case NPY_FLOAT:
          {
            float *values = reinterpret_cast<float*>(array->data);
            new (storage) Vector3x(values[0], values[1], values[2]);
          }
          break;
        case NPY_DOUBLE:
          {
            double *values = reinterpret_cast<double*>(array->data);
            new (storage) Vector3x(values[0], values[1], values[2]);
          }
          break;
        default:
          return;
      }

      data->convertible = storage;
    }
  };

void expose_eigen_vectors()
{
  import_array(); // needed for NumPy 

  
  Vector3x_to_python_array<math::rVector3d>();
  Vector3x_from_python_array<math::rVector3d>();
  Vector3x_to_python_array<math::iVector3d>();
  Vector3x_from_python_array<math::iVector3d>();
}
}}

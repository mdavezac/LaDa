// Adapted from Avogadro open source molecular editor from the kde project. (libavogadro/src/python/eigen.cpp)
#include <boost/python/detail/wrap_python.hpp>
#include <numpy/arrayobject.h> 
#include <boost/python.hpp>
#include <boost/python/tuple.hpp>

#include <Eigen/Geometry>

#include <python/std_vector.hpp>
#include "../eigen.h"

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
   * Matrix3x = Matrix3d, Matrix3f, Matrix3i
   *
   ***********************************************************************/
  
template <class Matrix3x>
  struct Matrix3x_to_python_array
  {
    typedef typename Matrix3x::Scalar Scalar;
    
    struct innerclass
    {
      //
      //  Eigen::Matrix3x --> python array
      //
      static PyObject* convert(Matrix3x const &vec)
      {
        npy_intp dims[2] = { 3, 3 };
        PyObject *result;
        if (ScalarTraits<Scalar>::isInt)
          result = PyArray_SimpleNew(2, dims, NPY_INT);
        else if (ScalarTraits<Scalar>::isFloat)
          result = PyArray_SimpleNew(2, dims, NPY_FLOAT);
        else
          result = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
        
        // copy the data
        Scalar *data = (Scalar*) reinterpret_cast<PyArrayObject*>(result)->data;
        for(size_t i(0); i < 3; ++i)
          for(size_t j(0); j < 3; ++j)
            data[i*3 + j] = vec(i,j);
        
        return incref(result);
      }
 
      //
      //  Eigen::Matrix3x * --> python array
      //    
      static PyObject* convert(Matrix3x *vec)
      {
        if (!vec) throw_error_already_set();
 
        npy_intp dims[2] = { 3, 3 };
        PyObject *result;
        if (ScalarTraits<Scalar>::isInt)
          result = PyArray_SimpleNew(2, dims, NPY_INT);
        else if (ScalarTraits<Scalar>::isFloat)
          result = PyArray_SimpleNew(2, dims, NPY_FLOAT);
        else
          result = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
        
        // copy the data
        Scalar *data = (Scalar*) reinterpret_cast<PyArrayObject*>(result)->data;
        for(size_t i(0); i < 3; ++i)
          for(size_t j(0); j < 3; ++j)
            data[i*3 + j] = (*vec)(i,j);
        
        return incref(result);
      }

      //
      //  const Eigen::Matrix3x * --> python array
      //
      static PyObject* convert(const Matrix3x *vec)
      {
        if (!vec)
          throw_error_already_set();

        npy_intp dims[2] = {3, 3};
        PyObject *result;
        if (ScalarTraits<Scalar>::isInt)
          result = PyArray_SimpleNew(2, dims, NPY_INT);
        else if (ScalarTraits<Scalar>::isFloat)
          result = PyArray_SimpleNew(2, dims, NPY_FLOAT);
        else
          result = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
        
        // copy the data
        Scalar *data = (Scalar*) reinterpret_cast<PyArrayObject*>(result)->data;
        for(size_t i(0); i < 3; ++i)
          for(size_t j(0); j < 3; ++j)
            data[i*3 + j] = (*vec)(i,j);
        
        return incref(result);
      }
    };
   
    Matrix3x_to_python_array()
    {
      to_python_converter<Matrix3x, innerclass>();
      to_python_converter<Matrix3x*, innerclass>();
      to_python_converter<const Matrix3x*, innerclass>();
    }

  };

  template <class Matrix3x>
  struct Matrix3x_from_python_array
  {
    typedef typename Matrix3x::Scalar Scalar;
    
    Matrix3x_from_python_array()
    {
      // Insert an rvalue from_python converter at the tail of the
      // chain. Used for implicit conversions
      //
      //  python array --> Matrix3x
      //
      // used for:
      //
      //  void function(Matrix3x vec)
      //  void function(Matrix3x & vec)
      //  void function(const Matrix3x & vec)
      //
      converter::registry::push_back( &convertible, &construct, type_id<Matrix3x>() );

      
      // Insert an lvalue from_python converter
      //
      //  python array --> Matrix3x*
      //
      // used for:
      //  
      //  void function(const Matrix3x * vec)
      converter::registry::insert( &convert, type_id<Matrix3x>() );
    }

    static void* convert(PyObject *obj_ptr)
    {
      if (!PyArray_Check(obj_ptr)) throw_error_already_set();

      // only accept int, long, float and double
      switch (PyArray_ObjectType(obj_ptr, 0)) 
      {
        case NPY_INT:
        case NPY_LONG:
        case NPY_FLOAT:
        case NPY_DOUBLE:
        break;
        default: PyErr_SetString(PyExc_ValueError, "Numpy array has wrong type.\n");
                 throw_error_already_set();
                 return NULL;
      }

      // do some type checking
      if ((PyArray_ObjectType(obj_ptr, 0) == NPY_FLOAT) || (PyArray_ObjectType(obj_ptr, 0) == NPY_DOUBLE))
        if (ScalarTraits<Scalar>::isInt)
        { 
          PyErr_SetString(PyExc_ValueError, "Numpy array has wrong type.\n");
          throw_error_already_set();
          return NULL;
        }

      if ((PyArray_ObjectType(obj_ptr, 0) == NPY_INT) || (PyArray_ObjectType(obj_ptr, 0) == NPY_LONG))
        if (ScalarTraits<Scalar>::isFloat || ScalarTraits<Scalar>::isDouble)
        { 
          PyErr_SetString(PyExc_ValueError, "Numpy array has wrong type.\n");
          throw_error_already_set();
          return NULL;
        }
      
      PyArrayObject *array = reinterpret_cast<PyArrayObject*>(obj_ptr);

      // check the dimensions
      if( array->nd == 1 and array->dimensions[0] != 9 )
      {
        PyErr_SetString(PyExc_ValueError, "1-dimensional numpy array cannot be converted to 3x3 Matrix.");
        throw_error_already_set(); // the array has at least two dimensions (matrix)
        return NULL;
      }
      else if (array->nd == 2 and (array->dimensions[0] != 3 or array->dimensions[1] != 3) )
      {
        PyErr_SetString(PyExc_ValueError, "2-dimensional numpy array cannot be converted to 3x3 Matrix.");
        throw_error_already_set(); // the array has at least two dimensions (matrix)
        return NULL;
      }
  

      switch (PyArray_ObjectType(obj_ptr, 0)) 
      {
        case NPY_INT: return create_matrix<int>( array, new Matrix3x );
        case NPY_LONG: return create_matrix<long>( array, new Matrix3x );
        case NPY_FLOAT: return create_matrix<float>( array, new Matrix3x );
        case NPY_DOUBLE: return create_matrix<double>( array, new Matrix3x );
        default: break;
      }
      return NULL;
    }

 
    static void* convertible(PyObject *obj_ptr)
    {
      if (!PyArray_Check(obj_ptr)) return 0;

      // only accept int, long, float and double
      switch (PyArray_ObjectType(obj_ptr, 0))
      {
        case NPY_INT:
        case NPY_LONG:
        case NPY_FLOAT:
        case NPY_DOUBLE: break;
        default: return NULL;
      }
      
      // do some type checking
      if ((PyArray_ObjectType(obj_ptr, 0) == NPY_FLOAT) || (PyArray_ObjectType(obj_ptr, 0) == NPY_DOUBLE))
        if (ScalarTraits<Scalar>::isInt) return 0;

      if ((PyArray_ObjectType(obj_ptr, 0) == NPY_INT) || (PyArray_ObjectType(obj_ptr, 0) == NPY_LONG))
        if (ScalarTraits<Scalar>::isFloat || ScalarTraits<Scalar>::isDouble) return 0;
      
      PyArrayObject *array = reinterpret_cast<PyArrayObject*>(obj_ptr);

      // check the dimensions
      if(    (array->nd == 1 and array->dimensions[0] != 3) 
          or (array->nd == 2 and (array->dimensions[0] != 3 or array->dimensions[1] != 3)) )
        return 0; 
      return obj_ptr;
    }

    static void construct(PyObject *obj_ptr, converter::rvalue_from_python_stage1_data *data)
    {
      PyArrayObject *array = reinterpret_cast<PyArrayObject*>(obj_ptr);
      void *storage = ((converter::rvalue_from_python_storage<Matrix3x>*)data)->storage.bytes;

      switch (PyArray_ObjectType(obj_ptr, 0)) 
      {
        case NPY_INT: create_matrix<int>( array, new (storage) Matrix3x ); break;
        case NPY_LONG: create_matrix<long>( array, new (storage) Matrix3x ); break;
        case NPY_FLOAT: create_matrix<float>( array, new (storage) Matrix3x ); break;
        case NPY_DOUBLE: create_matrix<double>( array, new (storage) Matrix3x ); break;
        default: return;
      }
      data->convertible = storage;
    }
    template<typename T2> 
      static void* create_matrix( PyArrayObject* array, Matrix3x *result )
      {
        T2 *values = reinterpret_cast<T2*>(array->data);
        if(array->nd == 1)
        {
          size_t strides = array->strides[0] / sizeof(T2);
          for(size_t i(0); i < 3; ++i)
            for(size_t j(0); j < 3; ++j)
              (*result)(i, j) = (Scalar)*(values + strides*(3*i+j));
        }
        else
        {
          size_t strides[2] = { array->strides[0] / sizeof(T2), array->strides[1] / sizeof(T2) };
          for(size_t i(0); i < 3; ++i)
            for(size_t j(0); j < 3; ++j)
              (*result)(i,j) = (Scalar)*(values + strides[0]*i + strides[1]*j);
        }
        return (void*)result;
      }

  };
  

void expose_eigen_matrices()
{
  import_array(); // needed for NumPy 

  
  Matrix3x_to_python_array<math::rMatrix3d>();
  Matrix3x_from_python_array<math::rMatrix3d>();
  Matrix3x_to_python_array<math::iMatrix3d>();
  Matrix3x_from_python_array<math::iMatrix3d>();

  expose_vector<math::rMatrix3d>("_rMatrix3dVector", "Exposes std::vector<math::rMatrix3d>.");
  expose_vector<math::iMatrix3d>("_iMatrix3dVector", "Exposes std::vector<math::rMatrix3d>.");
}
}}

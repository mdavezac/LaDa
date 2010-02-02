// Last update: timvdm 19 June 2009
// Adapted from Avogadro open source molecular editor from the kde project.
// (libavogadro/src/python/eigen.cpp)
#include <iostream>

#include <boost/type_traits/is_integral.hpp> 
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/python/detail/wrap_python.hpp>
#include <numpy/arrayobject.h> 
#include <boost/python.hpp>
#include <boost/python/tuple.hpp>

#include <Eigen/Geometry>

#include <python/std_vector.hpp>
#include "../eigen.h"
#include "numpy_types.h"
#include "avogadro.hpp"


namespace LaDa
{ 
  namespace python
  {
    using namespace boost::python;


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
          PyObject *result = PyArray_SimpleNew(1, dims, numpy::type<Scalar>::value);
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
          PyObject *result = PyArray_SimpleNew(1, dims, numpy::type<Scalar>::value);
          
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
          PyObject *result = PyArray_SimpleNew(1, dims, numpy::type<Scalar>::value);
          
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
        if (!PyArray_Check(obj_ptr)) throw_error_already_set();
  
        // only accept some type.
        if(    (boost::is_integral<Scalar>::value and not numpy::is_integer(obj_ptr))
            or (boost::is_floating_point<Scalar>::value and not numpy::is_float(obj_ptr)) )
        { 
          PyErr_SetString(PyExc_ValueError, "Numpy array has wrong type.\n");
          throw_error_already_set();
          return NULL;
        }
        
        PyArrayObject *array = reinterpret_cast<PyArrayObject*>(obj_ptr);
  
        // check the dimensions
        if (array->nd == 1 and array->dimensions[0] != 3)
        {
          PyErr_SetString(PyExc_ValueError, "Numpy array is 1d, but longer than 3.\n");
          throw_error_already_set(); // the array has at least two dimensions (matrix)
          return NULL;
        }
        else if(     array->nd ==  2
                 and (    (array->dimensions[0] == 1 and array->dimensions[1] != 3)
                       or (array->dimensions[0] == 3 and array->dimensions[1] != 1) ) )
        {
          PyErr_SetString(PyExc_ValueError, "Numpy array should be 3, 3x1, or 1x3.\n");
          throw_error_already_set(); // the 1D array does not have exactly 3 elements
          return NULL;
        }
  
        switch (PyArray_ObjectType(obj_ptr, 0)) 
        {
          case NPY_INT: return create_vector<int>( array, new Vector3x );
          case NPY_LONG: return create_vector<long>( array, new Vector3x );
          case NPY_FLOAT: return create_vector<float>( array, new Vector3x );
          case NPY_DOUBLE: return create_vector<double>( array, new Vector3x );
          default: break;
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
        if( boost::is_integral<Scalar>::value and not numpy::is_integer(obj_ptr) ) return NULL;
        if( boost::is_floating_point<Scalar>::value and not numpy::is_float(obj_ptr) ) return NULL;
        
        PyArrayObject *array = reinterpret_cast<PyArrayObject*>(obj_ptr);
  
        // check the dimensions
        if(array->nd < 1) return NULL;
        else if(array->nd == 1) { if( array->dimensions[0] != 3 ) return NULL; }
        else if(array->nd == 2)
        { 
          if( array->dimensions[0] == 1 and array->dimensions[1] != 3 ) return NULL; 
          else if( array->dimensions[0] == 3 and array->dimensions[1] != 1 ) return NULL; 
        }
        return obj_ptr;
      }
  
      static void construct(PyObject *obj_ptr, converter::rvalue_from_python_stage1_data *data)
      {
        PyArrayObject *array = reinterpret_cast<PyArrayObject*>(obj_ptr);
        void *storage = ((converter::rvalue_from_python_storage<Vector3x>*)data)->storage.bytes;
  
        switch (PyArray_ObjectType(obj_ptr, 0)) 
        {
          case NPY_INT: create_vector<int>( array, new (storage) Vector3x ); break;
          case NPY_LONG: create_vector<long>( array, new (storage) Vector3x ); break;
          case NPY_FLOAT: create_vector<float>( array, new (storage) Vector3x ); break;
          case NPY_DOUBLE: create_vector<double>( array, new (storage) Vector3x ); break;
          default: return;
        }
        data->convertible = storage;
      }
  
      template<typename T2> 
        static void* create_vector( PyArrayObject* array, Vector3x *result )
        {
          T2 *values = reinterpret_cast<T2*>(array->data);
          size_t const strides = array->strides[array->dimensions[0] == 3 ? 0: 1] / sizeof(T2);
#         ifdef _LADADEBUG
            if( array->strides[array->dimensions[0] == 3 ? 0: 1] % sizeof(T2) != 0 )
            {
              PyErr_SetString(PyExc_RuntimeError, "Incoherent numpy array.\n");
              bp::throw_error_already_set();
              return NULL;
            }
#         endif
          (*result)(0) = (Scalar) (*values);
          (*result)(1) = (Scalar) ( *(values + strides) );
          (*result)(2) = (Scalar) ( *(values + 2*strides) );
          return (void*)result;
        }
    };

    void expose_eigen_vectors()
    {
      import_array(); // needed for NumPy 
      
      Vector3x_to_python_array<math::rVector3d>();
      Vector3x_from_python_array<math::rVector3d>();
      Vector3x_to_python_array<math::iVector3d>();
      Vector3x_from_python_array<math::iVector3d>();
  
      expose_vector<math::rVector3d>("_rVector3dVector", "Exposes std::vector<math::rVector3d>.");
      expose_vector<math::iVector3d>("_iVector3dVector", "Exposes std::vector<math::iVector3d>.");
    }
  }
}

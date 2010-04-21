// Adapted from Avogadro open source molecular editor from the kde project. (libavogadro/src/python/eigen.cpp)
#include <iostream>

#include <boost/type_traits/is_integral.hpp> 
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/python/detail/wrap_python.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/errors.hpp>

#include <Eigen/Geometry>

#include <python/std_vector.hpp>
#include <python/numpy_types.h>

#include "../eigen.h"


namespace LaDa
{ 
  namespace python
  {
    using namespace boost::python;
    namespace numpy = LaDa::math::numpy;

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
            PyObject *result = PyArray_SimpleNew(2, dims, numpy::type<Scalar>::value);
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
            PyObject *result = PyArray_SimpleNew(2, dims, numpy::type<Scalar>::value);
            
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
            PyObject *result = PyArray_SimpleNew(2, dims, numpy::type<Scalar>::value);
            
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
          if(    (boost::is_integral<Scalar>::value and not numpy::is_integer(obj_ptr))
              or (boost::is_floating_point<Scalar>::value and not numpy::is_float(obj_ptr)) )
          { 
            PyErr_SetString(PyExc_ValueError, "Numpy array has wrong type.\n");
            throw_error_already_set();
            return NULL;
          }

          PyArrayObject *array = reinterpret_cast<PyArrayObject*>(obj_ptr);

          // check the dimensions
          if( array->nd == 1 and array->dimensions[0] != 9 )
          {
            PyErr_SetString(PyExc_ValueError,
                            "1-dimensional numpy array cannot be converted to 3x3 Matrix.");
            throw_error_already_set(); // the array has at least two dimensions (matrix)
            return NULL;
          }
          else if (array->nd == 2 and (array->dimensions[0] != 3 or array->dimensions[1] != 3) )
          {
            PyErr_SetString(PyExc_ValueError,
                            "2-dimensional numpy array cannot be converted to 3x3 Matrix.");
            throw_error_already_set(); // the array has at least two dimensions (matrix)
            return NULL;
          }
      

          switch (PyArray_ObjectType(obj_ptr, 0)) 
          {
            case NPY_FLOAT     : return create_matrix<npy_float     >(array, new Matrix3x);
            case NPY_DOUBLE    : return create_matrix<npy_double    >(array, new Matrix3x);
            case NPY_LONGDOUBLE: return create_matrix<npy_longdouble>(array, new Matrix3x);
            case NPY_INT       : return create_matrix<npy_int       >(array, new Matrix3x);
            case NPY_UINT      : return create_matrix<npy_uint      >(array, new Matrix3x);
            case NPY_LONG      : return create_matrix<npy_long      >(array, new Matrix3x);
            case NPY_ULONG     : return create_matrix<npy_ulong     >(array, new Matrix3x);
            case NPY_LONGLONG  : return create_matrix<npy_longlong  >(array, new Matrix3x);
            case NPY_ULONGLONG : return create_matrix<npy_ulonglong >(array, new Matrix3x);
            case NPY_BYTE      : return create_matrix<npy_byte      >(array, new Matrix3x);
            case NPY_UBYTE     : return create_matrix<npy_ubyte     >(array, new Matrix3x);
            case NPY_SHORT     : return create_matrix<npy_short     >(array, new Matrix3x);
            case NPY_USHORT    : return create_matrix<npy_ushort    >(array, new Matrix3x);
            default: break;
          }
          return NULL;
        }

     
        static void* convertible(PyObject *obj_ptr)
        {
          if (!PyArray_Check(obj_ptr)) return 0;

          // only accept int, long, float and double
          if(    (boost::is_integral<Scalar>::value and not numpy::is_integer(obj_ptr))
              or (boost::is_floating_point<Scalar>::value and not numpy::is_float(obj_ptr)) )
            return false;
          
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
            case NPY_FLOAT     : create_matrix<npy_float     >(array, new (storage) Matrix3x); break;
            case NPY_DOUBLE    : create_matrix<npy_double    >(array, new (storage) Matrix3x); break;
            case NPY_LONGDOUBLE: create_matrix<npy_longdouble>(array, new (storage) Matrix3x); break;
            case NPY_INT       : create_matrix<npy_int       >(array, new (storage) Matrix3x); break;
            case NPY_UINT      : create_matrix<npy_uint      >(array, new (storage) Matrix3x); break;
            case NPY_LONG      : create_matrix<npy_long      >(array, new (storage) Matrix3x); break;
            case NPY_ULONG     : create_matrix<npy_ulong     >(array, new (storage) Matrix3x); break;
            case NPY_LONGLONG  : create_matrix<npy_longlong  >(array, new (storage) Matrix3x); break;
            case NPY_ULONGLONG : create_matrix<npy_ulonglong >(array, new (storage) Matrix3x); break;
            case NPY_BYTE      : create_matrix<npy_byte      >(array, new (storage) Matrix3x); break;
            case NPY_UBYTE     : create_matrix<npy_ubyte     >(array, new (storage) Matrix3x); break;
            case NPY_SHORT     : create_matrix<npy_short     >(array, new (storage) Matrix3x); break;
            case NPY_USHORT    : create_matrix<npy_ushort    >(array, new (storage) Matrix3x); break;
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
    #         ifdef _LADADEBUG
                if( array->strides[0] % sizeof(T2) != 0 )
                {
                  PyErr_SetString(PyExc_RuntimeError, "Incoherent numpy array.\n");
                  throw_error_already_set();
                  return NULL;
                }
    #         endif
              for(size_t i(0); i < 3; ++i)
                for(size_t j(0); j < 3; ++j)
                  (*result)(i, j) = (Scalar)*(values + strides*(3*i+j));
            }
            else
            {
              size_t strides[2] = { array->strides[0] / sizeof(T2), array->strides[1] / sizeof(T2) };
    #         ifdef _LADADEBUG
                if( array->strides[0] % sizeof(T2) != 0 or array->strides[1] % sizeof(T2) )
                {
                  PyErr_SetString(PyExc_RuntimeError, "Incoherent numpy array.\n");
                  throw_error_already_set();
                  return NULL;
                }
    #         endif
              for(size_t i(0); i < 3; ++i)
                for(size_t j(0); j < 3; ++j)
                  (*result)(i,j) = (Scalar)*(values + strides[0]*i + strides[1]*j);
            }
            return (void*)result;
          }

      };
      

    void expose_eigen_matrices()
    {
      using namespace LaDa::Python;
      import_array(); // needed for NumPy 

      
      Matrix3x_to_python_array<math::rMatrix3d>();
      Matrix3x_from_python_array<math::rMatrix3d>();
      Matrix3x_to_python_array<math::iMatrix3d>();
      Matrix3x_from_python_array<math::iMatrix3d>();

      expose_vector<math::rMatrix3d>("_rMatrix3dVector", "Exposes std::vector<math::rMatrix3d>.");
      expose_vector<math::iMatrix3d>("_iMatrix3dVector", "Exposes std::vector<math::rMatrix3d>.");
    }
  }
}

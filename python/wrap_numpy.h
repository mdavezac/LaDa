#ifndef LADA_PYTHON_WRAP_NUMPY_H 
#define LADA_PYTHON_WRAP_NUMPY_H 

// Include numpy before including this file. 

#include "LaDaConfig.h"

#include <boost/type_traits/is_floating_point.hpp>

#include <python/object.h>
#include <math/eigen.h>
#include "numpy_types.h"

namespace LaDa
{
  namespace python
  {
    //! Convert/wrap a matrix to numpy.
    template<class T_DERIVED>
      PyObject* wrap_to_numpy(Eigen::DenseBase<T_DERIVED> const &_in, PyObject *_parent = NULL)
      {
        npy_intp dims[2] = { _in.rows(), _in.cols() };
        typedef math::numpy::type<typename Eigen::DenseBase<T_DERIVED>::Scalar> t_ScalarType;
        PyArrayObject *result = _parent == NULL ?
          (PyArrayObject*) PyArray_ZEROS(_in.cols() > 1? 2: 1, dims, t_ScalarType::value, _in.IsRowMajor?0:1):
          (PyArrayObject*) PyArray_SimpleNewFromData(_in.cols() > 1? 2: 1, dims, t_ScalarType::value,
                                     (void*)(&_in(0,0)));
        if(result == NULL) return NULL;
        // macro for row vs column major. The macro changed for npy >= 1.6
#       ifdef LADA_MACRO
#         error LADA_MACRO already defined
#       endif
#       ifdef NPY_ARRAY_C_CONTIGUOUS
#         define LADA_MACRO NPY_ARRAY_C_CONTIGUOUS;
#       else 
#         define LADA_MACRO NPY_C_CONTIGUOUS
#       endif
        // If has a parent, do not copy data, just incref it as base.
        if(_parent != NULL) 
        {
          // For some reason, eigen is column major, whereas c++ is generally row major.
          if(result->flags & LADA_MACRO and not _in.IsRowMajor) 
            result->flags -= LADA_MACRO;
          else if((not (result->flags & LADA_MACRO)) and _in.IsRowMajor) 
            result->flags |= LADA_MACRO;
          Eigen::DenseCoeffsBase<T_DERIVED> const coeffs = _in;
          if(_in.cols() == 1)
            result->strides[0] = _in.innerStride() * sizeof(typename t_ScalarType::np_type);
          else if(_in.IsRowMajor) 
          {
            result->strides[0] = _in.outerStride() * sizeof(typename t_ScalarType::np_type);
            result->strides[1] = _in.innerStride() * sizeof(typename t_ScalarType::np_type);
          }
          else 
          {
            result->strides[0] = _in.innerStride() * sizeof(typename t_ScalarType::np_type);
            result->strides[1] = _in.outerStride() * sizeof(typename t_ScalarType::np_type);
          }
          result->base = _parent;
          Py_INCREF(_parent);
        }
        // otherwise, copy data.
        else
        {
          for(size_t i(0); i < _in.rows(); ++i)
            for(size_t j(0); j < _in.cols(); ++j)
              *((typename t_ScalarType::np_type*)
                  (result->data + i*result->strides[0] + j*result->strides[1])) = _in(i, j);
        }
#       undef LADA_MACRO
#       ifdef NPY_ARRAY_WRITEABLE
#         define LADA_MACRO NPY_ARRAY_WRITEABLE
#       else
#         define LADA_MACRO NPY_WRITEABLE
#       endif
          if(result->flags & LADA_MACRO) result->flags -= LADA_MACRO;
#       undef LADA_MACRO
        return (PyObject*)result;
      }
    //! Convert/wrap a matrix to numpy.
    template<class T_DERIVED>
      PyObject* wrap_to_numpy(Eigen::DenseBase<T_DERIVED> &_in, PyObject *_parent = NULL)
      {
        npy_intp dims[2] = { _in.rows(), _in.cols() };
        typedef math::numpy::type<typename Eigen::DenseBase<T_DERIVED>::Scalar> t_ScalarType;
        PyArrayObject *result = _parent == NULL ?
          (PyArrayObject*) PyArray_ZEROS(_in.cols() > 1? 2: 1, dims, t_ScalarType::value, _in.IsRowMajor?0:1):
          (PyArrayObject*) PyArray_SimpleNewFromData(_in.cols() > 1? 2: 1, dims, t_ScalarType::value, &_in(0,0));
        if(result == NULL) return NULL;
        // macro for row vs column major. The macro changed for npy >= 1.6
#       ifdef LADA_MACRO
#         error LADA_MACRO already defined
#       endif
#       ifdef NPY_ARRAY_C_CONTIGUOUS
#         define LADA_MACRO NPY_ARRAY_C_CONTIGUOUS;
#       else 
#         define LADA_MACRO NPY_C_CONTIGUOUS
#       endif
        // If has a parent, do not copy data, just incref it as base.
        if(_parent != NULL) 
        {
          // For some reason, eigen is column major, whereas c++ is generally row major.
          if(result->flags & LADA_MACRO and not _in.IsRowMajor) 
            result->flags -= LADA_MACRO;
          else if((not (result->flags & LADA_MACRO)) and _in.IsRowMajor) 
            result->flags |= LADA_MACRO;
          Eigen::DenseCoeffsBase<T_DERIVED> const coeffs = _in;
          if(_in.cols() == 1)
            result->strides[0] = _in.innerStride() * sizeof(typename t_ScalarType::np_type);
          else if(_in.IsRowMajor) 
          {
            result->strides[0] = _in.outerStride() * sizeof(typename t_ScalarType::np_type);
            result->strides[1] = _in.innerStride() * sizeof(typename t_ScalarType::np_type);
          }
          else 
          {
            result->strides[0] = _in.innerStride() * sizeof(typename t_ScalarType::np_type);
            result->strides[1] = _in.outerStride() * sizeof(typename t_ScalarType::np_type);
          }
          result->base = _parent;
          Py_INCREF(_parent);
        }
        // otherwise, copy data.
        else
        {
          for(size_t i(0); i < _in.rows(); ++i)
            for(size_t j(0); j < _in.cols(); ++j)
              *((typename t_ScalarType::np_type*)
                  (result->data + i*result->strides[0] + j*result->strides[1])) = _in(i, j);
        }
#       undef LADA_MACRO
        return (PyObject*)result;
      }
    //! Converts an input sequence to a cell.
    template<class T_DERIVED>
      bool convert_to_matrix(PyObject *_in, Eigen::DenseBase<T_DERIVED> &_out)
      {
        size_t const N0(_out.rows());
        size_t const N1(_out.cols());
        if(PyArray_Check(_in))
        {
          if(PyArray_NDIM(_in) != 2)
          {
            npy_intp const n(PyArray_NDIM(_in));
            LADA_PYERROR_FORMAT(TypeError, "Expected a 2d array, got %id", n);
            return false;
          }
          if(PyArray_DIM(_in, 0) != N0 or PyArray_DIM(_in, 1) != N1)
          {
            npy_intp const n0(PyArray_DIM(_in, 0));
            npy_intp const n1(PyArray_DIM(_in, 1));
            PyErr_Format( ::LaDa::python::PyException< ::LaDa::error::TypeError >::exception().ptr(),
                          "Expected a %ix%i array, got %ix%i", N0, N1, n0, n1);
            return false;
          }
          python::Object iterator = PyArray_IterNew(_in);
          if(not iterator) return false;
          int const type = PyArray_DESCR(_in)->type_num;
#         ifdef  LADA_NPYITER
#           error LADA_NPYITER already defined
#         endif
#         define LADA_NPYITER(TYPE, NUM_TYPE)                                             \
            if(type == NUM_TYPE)                                                          \
            {                                                                             \
              for(size_t i(0); i < N0*N1; ++i)                                            \
              {                                                                           \
                if(not PyArray_ITER_NOTDONE(iterator.borrowed()))                         \
                {                                                                         \
                  LADA_PYERROR(TypeError, "Numpy array too small.");                      \
                  return false;                                                           \
                }                                                                         \
                _out(i/N1, i%N1) = *((TYPE*) PyArray_ITER_DATA(iterator.borrowed()));     \
                PyArray_ITER_NEXT(iterator.borrowed());                                   \
              }                                                                           \
              if(PyArray_ITER_NOTDONE(iterator.borrowed()))                               \
              {                                                                           \
                LADA_PYERROR(TypeError, "Numpy array too long.");                         \
                return false;                                                             \
              }                                                                           \
            }
          LADA_NPYITER( npy_float,      NPY_FLOAT)      
          else LADA_NPYITER( npy_double,     NPY_DOUBLE     )
          else LADA_NPYITER( npy_longdouble, NPY_LONGDOUBLE )
          else LADA_NPYITER( npy_int,        NPY_INT        )
          else LADA_NPYITER( npy_uint,       NPY_UINT       )
          else LADA_NPYITER( npy_long,       NPY_LONG       )
          else LADA_NPYITER( npy_longlong,   NPY_LONGLONG   )
          else LADA_NPYITER( npy_ulonglong,  NPY_ULONGLONG  )
          else LADA_NPYITER( npy_ubyte,      NPY_BYTE       )
          else LADA_NPYITER( npy_short,      NPY_SHORT      )
          else LADA_NPYITER( npy_ushort,     NPY_USHORT     )
          else
          {
            LADA_PYERROR(TypeError, "Unknown numpy array type.");
            return false;
          }
#         undef LADA_NPYITER
        } // numpy array
        else 
        {
          python::Object i_outer = PyObject_GetIter(_in);
          if(not i_outer) { return false; }
          python::Object outer(PyIter_Next(i_outer.borrowed()));
          if(not outer.is_valid()) return false; 
          if(not outer.hasattr("__iter__")) // except 9 in a row
          {
            size_t i(0);
            for( ; outer.is_valid() and i < N0*N1;
                 outer.reset(PyIter_Next(i_outer.borrowed())), ++i ) 
            {
              if(PyInt_Check(outer.borrowed())) _out(i/N0, i%N0) = PyInt_AS_LONG(outer.borrowed());
              else if(PyFloat_Check(outer.borrowed())) _out(i/N0, i%N0) = PyFloat_AS_DOUBLE(outer.borrowed());
              else
              { 
                LADA_PYERROR(TypeError, "Object should contains numbers only.");
                return false;
              }
            }
            if(outer.is_valid() or i != 9)
            {
              LADA_PYERROR(TypeError, "Expected 9 (NxN) numbers in input.");
              return false;
            }
          }    // N0*N1 in a row.
          else // expect N0 by N1
          {
            size_t i(0);
            for( ; outer.is_valid() and i < N0;
                 outer.reset(PyIter_Next(i_outer.borrowed())), ++i ) 
            {
              python::Object i_inner = PyObject_GetIter(outer.borrowed());
              if(not i_inner) return false;
              python::Object inner(PyIter_Next(i_inner.borrowed()));
              if(not inner) return false;
              size_t j(0);
              for( ; inner.is_valid() and j < N1;
                   inner.reset(PyIter_Next(i_inner.borrowed())), ++j ) 
              {
                if(PyInt_Check(inner.borrowed())) _out(i, j) = PyInt_AS_LONG(inner.borrowed());
                else if(PyFloat_Check(inner.borrowed())) _out(i, j) = PyFloat_AS_DOUBLE(inner.borrowed());
                else
                { 
                  LADA_PYERROR(TypeError, "Object should contains numbers only.");
                  return false;
                }
              } // inner loop.
              if(inner.is_valid() or j != N1)
              {
                LADA_PYERROR(TypeError, "Not a NxN matrix of numbers.");
                return false;
              }
            } // outer loop.
            if(outer.is_valid() or i != N1)
            {
              LADA_PYERROR(TypeError, "Not a NxN matrix of numbers.");
              return false;
            }
          }
        } // sequence.
        return true;
      }
    //! Converts an input sequence to a cell.
    template<class T_DERIVED> 
      bool convert_to_vector(PyObject *_in, Eigen::DenseBase<T_DERIVED> &_out)
      {
        if(PyArray_Check(_in))
        {
          python::Object iterator = PyArray_IterNew(_in);
          if(not iterator) return false;
          int const type = PyArray_DESCR(_in)->type_num;
#         ifdef LADA_NPYITER
#           error LADA_NPYITER is already defined.
#         endif
#         define LADA_NPYITER(TYPE, NUM_TYPE)                                        \
            if(type == NUM_TYPE)                                                     \
            {                                                                        \
              for(size_t i(0); i < 3; ++i)                                           \
              {                                                                      \
                if(not PyArray_ITER_NOTDONE(iterator.borrowed()))                    \
                {                                                                    \
                  LADA_PYERROR(TypeError, "Numpy array too small.");                 \
                  return false;                                                      \
                }                                                                    \
                _out[i] = *((TYPE*) PyArray_ITER_DATA(iterator.borrowed()));         \
                PyArray_ITER_NEXT(iterator.borrowed());                              \
              }                                                                      \
              if(PyArray_ITER_NOTDONE(iterator.borrowed()))                          \
              {                                                                      \
                LADA_PYERROR(TypeError, "Numpy array too long.");                    \
                return false;                                                        \
              }                                                                      \
            }
          LADA_NPYITER( npy_float,      NPY_FLOAT)      
          else LADA_NPYITER( npy_double,     NPY_DOUBLE     )
          else LADA_NPYITER( npy_longdouble, NPY_LONGDOUBLE )
          else LADA_NPYITER( npy_int,        NPY_INT        )
          else LADA_NPYITER( npy_uint,       NPY_UINT       )
          else LADA_NPYITER( npy_long,       NPY_LONG       )
          else LADA_NPYITER( npy_longlong,   NPY_LONGLONG   )
          else LADA_NPYITER( npy_ulonglong,  NPY_ULONGLONG  )
          else LADA_NPYITER( npy_ubyte,      NPY_BYTE       )
          else LADA_NPYITER( npy_short,      NPY_SHORT      )
          else LADA_NPYITER( npy_ushort,     NPY_USHORT     )
          else
          {
            LADA_PYERROR(TypeError, "Unknown numpy array type.");
            return false;
          }
      #   undef LADA_NPYITER
        }
        else if(PyInt_Check(_in)) _out = Eigen::DenseBase<T_DERIVED>::Ones() * PyInt_AS_LONG(_in); 
        else if( boost::is_floating_point< typename Eigen::DenseBase<T_DERIVED> >::value 
                 and PyFloat_Check(_in))
          _out = Eigen::DenseBase<T_DERIVED>::Ones() * PyFloat_AS_DOUBLE(_in); 
        else
        {
          python::Object i_outer = PyObject_GetIter(_in);
          if(not i_outer) { return false; }
          python::Object item(PyIter_Next(i_outer.borrowed()));
          size_t i(0);
          for( ; item.is_valid() and i < 3;
               item.reset(PyIter_Next(i_outer.borrowed())), ++i ) 
            if(PyInt_Check(item.borrowed()) == 1) _out[i] = PyInt_AS_LONG(item.borrowed());
            else if(PyFloat_Check(item.borrowed()) == 1) _out[i] = PyFloat_AS_DOUBLE(item.borrowed());
            else
            { 
              LADA_PYERROR(TypeError, "Input vector should contains numbers only.");
              return false;
            }
          if(item.is_valid())
          { 
            LADA_PYERROR(TypeError, "Input vector is too large.");
            return false;
          }
          else if (i != 3) 
          {
            LADA_PYERROR(TypeError, "Input vector is too small.");
            return false;
          }
        } 
        return true;
      }
  }
}

#endif

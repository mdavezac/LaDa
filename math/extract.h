#ifndef LADA_MATH_EXTRACT 
#define LADA_MATH_EXTRACT 

// Include numpy before including this file. 

#include "LaDaConfig.h"

#include <python/object.h>

#include "eigen.h"

namespace LaDa
{
  namespace math
  {
    //! Converts an input sequence to a cell.
    template<class T_DERIVED>
      bool convert_to_cell(PyObject *_in, Eigen::DenseBase<T_DERIVED> &_out)
      {
        if(PyArray_Check(_in))
        {
          python::Object iterator = PyArray_IterNew(_in);
          if(not iterator) return false;
          int const type = PyArray_DESCR(_in)->type_num;
#         ifdef LADA_NPYITER
#           error LADA_NPYITER is already defined.
#         endif
#         define LADA_NPYITER(TYPE, NUM_TYPE)                                             \
            if(type == NUM_TYPE)                                                          \
            {                                                                             \
              for(size_t i(0); i < 9; ++i)                                                \
              {                                                                           \
                if(not PyArray_ITER_NOTDONE(iterator.borrowed()))                         \
                {                                                                         \
                  LADA_PYERROR(TypeError, "Numpy array too small.");                      \
                  return false;                                                           \
                }                                                                         \
                _out(i/3, i%3) = *((TYPE*) PyArray_ITER_DATA(iterator.borrowed()));       \
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
            for( ; outer.is_valid() and i < 9;
                 outer.reset(PyIter_Next(i_outer.borrowed())), ++i ) 
            {
              if(PyInt_Check(outer.borrowed())) _out(i/3, i%3) = PyInt_AS_LONG(outer.borrowed());
              else if(PyFloat_Check(outer.borrowed())) _out(i/3, i%3) = PyFloat_AS_DOUBLE(outer.borrowed());
              else
              { 
                LADA_PYERROR(TypeError, "Object should contains numbers only.");
                return false;
              }
            }
            if(outer.is_valid() or i != 9)
            {
              LADA_PYERROR(TypeError, "Expected 9 (3x3) numbers in input.");
              return false;
            }
          }    // 9 in a row.
          else // expect 3 by 3
          {
            size_t i(0);
            for( ; outer.is_valid() and i < 3;
                 outer.reset(PyIter_Next(i_outer.borrowed())), ++i ) 
            {
              python::Object i_inner = PyObject_GetIter(outer.borrowed());
              if(not i_inner) return false;
              python::Object inner(PyIter_Next(i_inner.borrowed()));
              if(not inner) return false;
              size_t j(0);
              for( ; inner.is_valid() and j < 3;
                   inner.reset(PyIter_Next(i_inner.borrowed())), ++j ) 
              {
                if(PyInt_Check(inner.borrowed())) _out(j, i) = PyInt_AS_LONG(inner.borrowed());
                else if(PyFloat_Check(inner.borrowed())) _out(j, i) = PyFloat_AS_DOUBLE(inner.borrowed());
                else
                { 
                  LADA_PYERROR(TypeError, "Object should contains numbers only.");
                  return false;
                }
              } // inner loop.
              if(inner.is_valid() or j != 3)
              {
                LADA_PYERROR(TypeError, "Not a 3x3 matrix of numbers.");
                return false;
              }
            } // outer loop.
            if(outer.is_valid() or i != 3)
            {
              LADA_PYERROR(TypeError, "Not a 3x3 matrix of numbers.");
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
        else if(PyInt_Check(_in)) _out = math::rVector3d::Ones() * PyInt_AS_LONG(_in); 
        else if(PyFloat_Check(_in)) _out = math::rVector3d::Ones() * PyFloat_AS_DOUBLE(_in); 
        else
        {
          python::Object i_outer = PyObject_GetIter(_in);
          if(not i_outer) { return false; }
          python::Object item(PyIter_Next(i_outer.borrowed()));
          size_t i(0);
          for( ; item.is_valid() or i >= 3;
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

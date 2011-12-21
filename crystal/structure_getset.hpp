namespace LaDa
{
  namespace crystal
  {
    extern "C"
    {
      //! Returns cell as a numpy array. 
      static PyObject* structure_getcell(StructureData *_self, void *closure);
      //! Sets cell from a sequence of 3x3 numbers.
      static int structure_setcell(StructureData *_self, PyObject *_value, void *_closure);
      // Returns the scale.
      static PyObject* structure_getscale(StructureData *_self, void *closure);
      //! Sets the scale from a number.
      static int structure_setscale(StructureData *_self, PyObject *_value, void *_closure);
    }
  
    // Returns cell as a numpy array. 
    // Numpy does not implement python's cyclic garbage, hence new wrapper need be
    // created each call.
    static PyObject* structure_getcell(StructureData *_self, void *closure)
    {
      npy_intp dims[2] = {3, 3};
      int const value = math::numpy::type<math::rMatrix3d::Scalar>::value;
      PyArrayObject *result = (PyArrayObject*) PyArray_SimpleNewFromData(2, dims, value, _self->cell.data());
      if(result == NULL) return NULL;
      result->base = (PyObject*)_self;
      Py_INCREF(_self); // Increfed as base of array.
      return (PyObject*)result;
    }
    // Sets cell from a sequence of three numbers.
    static int structure_setcell(StructureData *_self, PyObject *_value, void *_closure)
    {
      if(_value == NULL)
      {
        LADA_PYERROR(TypeError, "Cannot delete cell attribute.");
        return -1;
      }
      if(PyArray_Check(_value))
      {
        PyObject* iterator = PyArray_IterNew(_value);
        if(iterator == NULL) return -1;
        int const type = PyArray_DESCR(_value)->type_num;
#       ifdef LADA_NPYITER
#         error LADA_NPYITER is already defined.
#       endif
#       define LADA_NPYITER(TYPE, NUM_TYPE)                                             \
          if(type == NUM_TYPE)                                                          \
          {                                                                             \
            for(size_t i(0); i < 9; ++i)                                                \
            {                                                                           \
              if(not PyArray_ITER_NOTDONE(iterator))                                    \
              {                                                                         \
                Py_DECREF(iterator);                                                    \
                LADA_PYERROR(TypeError, "Numpy array too small.");                      \
                return -1;                                                              \
              }                                                                         \
              _self->cell(i/3, i%3) = *((TYPE*) PyArray_ITER_DATA(iterator)); \
              PyArray_ITER_NEXT(iterator);                                              \
            }                                                                           \
            if(PyArray_ITER_NOTDONE(iterator))                                          \
            {                                                                           \
              Py_DECREF(iterator);                                                      \
              LADA_PYERROR(TypeError, "Numpy array too long.");                         \
              return -1;                                                                \
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
          Py_DECREF(iterator); 
          LADA_PYERROR(TypeError, "Unknown numpy array type.");
          return -1;
        }
#       undef LADA_NPYITER
        Py_DECREF(iterator); 
      }
      else 
      {
        PyObject* i_outer = PyObject_GetIter(_value);
        if(_value == NULL) { return -1; }
        size_t i(0);
        while(PyObject* outer = PyIter_Next(i_outer))
        {
          if(i >= 3) 
          {
            LADA_PYERROR(TypeError, "Not a 3x3 matrix of numbers.");
            goto outer_error;
          }
          if(PyObject *i_inner = PyObject_GetIter(outer))
          {
            size_t j(0);
            while(PyObject* inner = PyIter_Next(i_inner))
            {
              if(j >= 3) 
              {
                LADA_PYERROR(TypeError, "Not a 3x3 matrix of numbers.");
                goto inner_error;
              }
              if(PyInt_Check(inner)) _self->cell(j, i) = PyInt_AS_LONG(inner);
              else if(PyFloat_Check(inner)) _self->cell(j, i) = PyFloat_AS_DOUBLE(inner);
              else
              { 
                LADA_PYERROR(TypeError, "Object should contains numbers only.");
                goto inner_error;
              }
              Py_DECREF(inner);
              ++j;
              continue;
              inner_error:
                Py_DECREF(i_outer);
                Py_DECREF(outer);
                Py_DECREF(i_inner);
                Py_DECREF(inner);
                return -1;
            }
            Py_DECREF(i_inner);
          }
          else
          { 
            LADA_PYERROR(TypeError, "Object should contains numbers only.");
            goto outer_error;
          }
          Py_DECREF(outer);
          ++i;
          continue;
          outer_error:
            Py_DECREF(i_outer);
            Py_DECREF(outer);
            return -1;
        }
        Py_DECREF(i_outer);
      } 
      return 0;
    }
  
    // Returns the scale of the structure.
    static PyObject* structure_getscale(StructureData *_self, void *closure)
      { return PyFloat_FromDouble(_self->scale); }
    // Sets the scale of the structure from a number.
    static int structure_setscale(StructureData *_self, PyObject *_value, void *_closure)
    {
      if(_value == NULL) 
      {
        LADA_PYERROR(TypeError, "Cannot delete scale attribute.");
        return -1;
      }
      if(python::Quantity::isinstance(_value))
        try { _self->scale = python::Quantity(_value).get("angstrom"); }
        catch(...) { return -1; }
      else if(PyFloat_Check(_value)) _self->scale = PyFloat_AS_DOUBLE(_value); 
      else if(PyInt_Check(_value)) _self->scale = PyInt_AS_LONG(_value); 
      else
      {
        LADA_PYERROR(TypeError, "Input to scale is not an acceptable type.");
        return -1;
      }
      return 0;
    }
  }
} 

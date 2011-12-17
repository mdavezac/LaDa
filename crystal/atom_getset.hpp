namespace LaDa
{
  namespace crystal
  {
    extern "C"
    {
      //! Returns position as a numpy array. 
      static PyObject* lada_atom_getpos(AtomData *_self, void *closure);
      //! Sets position from a sequence of three numbers.
      static int lada_atom_setpos(AtomData *_self, PyObject *_value, void *_closure);
      //! Returns type python object.
      static PyObject* lada_atom_gettype(AtomData *_self, void *closure);
      //! Sets type python object.
      static int lada_atom_settype(AtomData *_self, PyObject *_value, void *_closure);
    }
  
    // Returns position as a numpy array. 
    // Numpy does not implement python's cyclic garbage, hence new wrapper need be
    // created each call.
    static PyObject* lada_atom_getpos(AtomData *_self, void *closure)
    {
      npy_intp dims[1] = {3};
      int const value = math::numpy::type<math::rVector3d::Scalar>::value;
      PyArrayObject* result = (PyArrayObject*) PyArray_SimpleNewFromData(1, dims, value, _self->pos.data());
      if(result == NULL) return NULL;
      result->base = (PyObject*)_self;
      Py_INCREF(_self); // Increfed as base of array.
      return (PyObject*)result;
    }
    // Sets position from a sequence of three numbers.
    static int lada_atom_setpos(AtomData *_self, PyObject *_value, void *_closure)
    {
      if(_value == NULL)
      {
        LADA_PYERROR(TypeError, "Cannot delete pos attribute.");
        return -1;
      }
      if(PyArray_Check(_value))
      {
        PyObject* iterator = PyArray_IterNew(_value);
        if(iterator == NULL) return -1;
        int const type = PyArray_DESCR(_value)->type_num;
    #   ifdef LADA_NPYITER
    #     error LADA_NPYITER is already defined.
    #   endif
    #   define LADA_NPYITER(TYPE, NUM_TYPE)                                        \
          if(type == NUM_TYPE)                                                     \
          {                                                                        \
            for(size_t i(0); i < 3; ++i)                                           \
            {                                                                      \
              if(not PyArray_ITER_NOTDONE(iterator))                               \
              {                                                                    \
                Py_DECREF(iterator);                                               \
                LADA_PYERROR(TypeError, "Numpy array too small.");                 \
                return -1;                                                         \
              }                                                                    \
              _self->pos[i] = *((TYPE*) PyArray_ITER_DATA(iterator));        \
              PyArray_ITER_NEXT(iterator);                                         \
            }                                                                      \
            if(PyArray_ITER_NOTDONE(iterator))                                     \
            {                                                                      \
              Py_DECREF(iterator);                                                 \
              LADA_PYERROR(TypeError, "Numpy array too long.");                    \
              return -1;                                                           \
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
          Py_DECREF(iterator); 
          LADA_PYERROR(TypeError, "Unknown numpy array type.");
          return -1;
        }
    #   undef LADA_NPYITER
        Py_DECREF(iterator); 
      }
      else if(PyInt_Check(_value)) _self->pos = math::rVector3d::Ones() * PyInt_AS_LONG(_value); 
      else if(PyFloat_Check(_value)) _self->pos = math::rVector3d::Ones() * PyFloat_AS_DOUBLE(_value); 
      else
      {
        PyObject* iterator = PyObject_GetIter(_value);
        if(iterator == NULL) { return -1; }
        size_t i(0);
        while(PyObject* item = PyIter_Next(iterator))
        {
          if(i >= 3)
          {
            LADA_PYERROR(TypeError, "Object is too large.");
            Py_DECREF(iterator);
            Py_DECREF(item);
            return -1; 
          } 
          if(PyInt_Check(item) == 1) _self->pos[i] = PyInt_AS_LONG(item);
          else if(PyFloat_Check(item) == 1) _self->pos[i] = PyFloat_AS_DOUBLE(item);
          else
          { 
            LADA_PYERROR(TypeError, "Object should contains numbers only.");
            Py_DECREF(iterator);
            Py_DECREF(item);
            return -1;
          }
          Py_DECREF(item);
          ++i;
        }
        if(i != 3)
        { 
          LADA_PYERROR(TypeError, "Object is too small.");
          Py_DECREF(iterator);
          return -1; 
        }
        Py_DECREF(iterator);
      } 
      return 0;
    }
    // Returns type python object.
    PyObject* lada_atom_gettype(AtomData *_self, void *closure)
      { Py_INCREF(_self->type); return _self->type; }
      
    // Sets type python object.
    int lada_atom_settype(AtomData *_self, PyObject *_value, void *_closure)
    {
      if(_value == NULL)
      {
        LADA_PYERROR(TypeError, "Cannot delete type. You can set it to None though.");
        return -1;
      }
      PyObject* dummy = _self->type;
      _self->type = _value;
      Py_INCREF(_value);
      Py_DECREF(dummy);
      return 0;
    }
  } // namespace Crystal
} // namespace LaDa
  

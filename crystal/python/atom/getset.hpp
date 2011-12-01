extern "C"
{
  //! Returns position as a numpy array. 
 static PyObject* LADA_NAME(getpos)(LADA_TYPE *_self, void *closure);
  //! Sets position from a sequence of three numbers.
  static int LADA_NAME(setpos)(LADA_TYPE *_self, PyObject *_value, void *_closure);
  //! Returns the atomic type. 
  static PyObject* LADA_NAME(getsite)(LADA_TYPE *_self, void *closure);
  //! Sets the atomic type.
  static int LADA_NAME(setsite)(LADA_TYPE *_self, PyObject *_value, void *_closure);
  // Returns the freeze flag.
  static PyObject* LADA_NAME(getfreeze)(LADA_TYPE *_self, void *closure);
  //! Sets the freeze flag from an unsigned.
  static int LADA_NAME(setfreeze)(LADA_TYPE *_self, PyObject *_value, void *_closure);
  // Returns the type of the atom.
  static PyObject* LADA_NAME(gettype)(LADA_TYPE *_self, void *closure);
  //! Sets the type from a string.
  static int LADA_NAME(settype)(LADA_TYPE *_self, PyObject *_value, void *_closure);
  //! Gets dictionary.
  static PyObject* LADA_NAME(getdict)(LADA_TYPE *_self, void *_closure);
  //! Sets dictionary.
  static int LADA_NAME(setdict)(LADA_TYPE *_self, PyObject *_value, void *_closure);
}

// Returns position as a numpy array. 
// Numpy does not implement python's cyclic garbage, hence new wrapper need be
// created each call.
static PyObject* LADA_NAME(getpos)(LADA_TYPE *_self, void *closure)
{
  npy_intp dims[1] = {3};
  int const value = math::numpy::type<math::rVector3d::Scalar>::value;
  PyArrayObject* result = (PyArrayObject*) PyArray_SimpleNewFromData(1, dims, value, _self->atom->pos.data());
  if(result == NULL) return NULL;
  result->base = (PyObject*)_self;
  Py_INCREF(_self); // Increfed as base of array.
  return (PyObject*)result;
}
// Sets position from a sequence of three numbers.
static int LADA_NAME(setpos)(LADA_TYPE *_self, PyObject *_value, void *_closure)
{
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
          _self->atom->pos[i] = *((TYPE*) PyArray_ITER_DATA(iterator));        \
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
  else if(PyInt_Check(_value)) _self->atom->pos = math::rVector3d::Ones() * PyInt_AS_LONG(_value); 
  else if(PyFloat_Check(_value)) _self->atom->pos = math::rVector3d::Ones() * PyFloat_AS_DOUBLE(_value); 
  else
  {
    PyObject* iterator = PyObject_GetIter(_value);
    if(_value == NULL) { return -1; }
    for(size_t i(0); i < 3; ++i)
      if(PyObject* item = PyIter_Next(iterator))
      {
        if(PyInt_Check(item) == 1) _self->atom->pos[i] = PyInt_AS_LONG(item);
        else if(PyFloat_Check(item) == 1) _self->atom->pos[i] = PyFloat_AS_DOUBLE(item);
        else
        { 
          LADA_PYERROR(TypeError, "Object should contains numbers only.");
          Py_DECREF(iterator);
          Py_DECREF(item);
          return -1;
        }
        Py_DECREF(item);
      }
      else
      { 
        LADA_PYERROR(TypeError, "Object is too small.");
        Py_DECREF(iterator);
        return -1; 
      }
    if(PyObject* item = PyIter_Next(iterator))
    {
      LADA_PYERROR(TypeError, "Object is too long.");
      Py_DECREF(item);
      Py_DECREF(iterator);
      return -1; 
    }
    Py_DECREF(iterator);
  } 
  return 0;
}

// Returns the atomic type. 
static PyObject* LADA_NAME(getsite)(LADA_TYPE *_self, void *closure)
  { return PyInt_FromLong(_self->atom->site); }
// Sets the atomic type.
static int LADA_NAME(setsite)(LADA_TYPE *_self, PyObject *_value, void *_closure)
{
  long const result = PyInt_AsLong(_value);
  if(result == -1 and PyErr_Occurred() != NULL) return -1;
  _self->atom->site = result;
  return 0;
}
// Returns the freeze flag.
static PyObject* LADA_NAME(getfreeze)(LADA_TYPE *_self, void *closure)
{
  long result = _self->atom->freeze;
  return PyInt_FromLong(result);
}
// Sets the freeze flag from an unsigned.
static int LADA_NAME(setfreeze)(LADA_TYPE *_self, PyObject *_value, void *_closure)
{
  long const result = PyInt_AsLong(_value);
  if(result == -1 and PyErr_Occurred() != NULL) return -1;
  if(result < 0)
  {
    PyErr_SetString( PyException<error::ValueError>::exception().ptr(),  
                     "Cannot set freeze to a negative value." );
    return -1;
  }
  _self->atom->freeze = result;
  return 0;
}
# if LADA_ATOM_NUMBER == 0
  static PyObject *LADA_NAME(gettype)(LADA_TYPE *_self, void *closure)
    { return PyString_FromString(_self->atom->type.c_str()); }
  static int LADA_NAME(settype)(LADA_TYPE *_self, PyObject *_value, void *closure)
  {
    if(char * const string = PyString_AsString(_value))
    {
      _self->atom->type = string;
      return 0;
    }
    return -1;
  }
# elif LADA_ATOM_NUMBER == 1
  static PyObject *LADA_NAME(gettype)(LADA_TYPE *_self, void *closure) 
    { Py_INCREF(_self->sequence); return (PyObject*) sequence; }
  static int LADA_NAME(settype)(LADA_TYPE *_self, PyObject *_value, void *closure)
    { to_cpp_sequence(_value, _self->atom->type); return 0;}
# endif
// Gets dictionary.
static PyObject* LADA_NAME(getdict)(LADA_TYPE *_self, void *_closure)
  { Py_INCREF(_self->dictionary); return _self->dictionary; }
// Sets dictionary.
static int LADA_NAME(setdict)(LADA_TYPE *_self, PyObject *_value, void *_closure)
{
  PyObject* dummy = _self->dictionary;
  _self->dictionary = _value;
  Py_DECREF(dummy);
  return 0;
}

static char posdoc[] = "Position in cartesian coordinates.";
static char posname[] = "pos";
static char sitedoc[] = "Site index (integer).\n\n"
                        "Generally given for a supercell with respect to a backbone lattice."
                        "It is -1 if not initialized, and positive or zero otherwise."
                        "It refers to the original atom in the lattice.";
static char dict_doc[] = "Atom dictionary holding python attributes.";
# if LADA_ATOM_NUMBER  == 0
  static char typedoc[] = "Atomic specie. Must be a string.";
# elif LADA_ATOM_NUMBER == 1
  static char typedoc[] = "List of atomic species.";
# endif 
static char sitename[] = "site";
static char freezedoc[] = "Mask to freeze position or type (unsigned integer).";
static char freezename[] = "freeze";
static char type_name[] = "type";
static char dict_name[] = "__dict__";

extern "C" 
{
  static PyGetSetDef LADA_NAME(getsetters)[] = {
      {posname, (getter)LADA_NAME(getpos), (setter)LADA_NAME(setpos), posdoc, NULL},
      {sitename, (getter)LADA_NAME(getsite), (setter)LADA_NAME(setsite), sitedoc, NULL},
      {freezename, (getter)LADA_NAME(getfreeze), (setter)LADA_NAME(setfreeze), freezedoc, NULL}, 
      {type_name, (getter)LADA_NAME(gettype), (setter)LADA_NAME(settype), typedoc, NULL}, 
      {dict_name, (getter)LADA_NAME(getdict), (setter)LADA_NAME(setdict), typedoc, NULL}, 
      {NULL}  /* Sentinel */
  };
}

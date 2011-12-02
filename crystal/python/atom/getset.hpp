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
  //! \brief Gets an attribute.
  //! \details Calls python's _PyObject_GenericGetAttrWithDict(...) with dict
  //!          set to the one contained in atom. This way, we are sure there is
  //!          a single dictionary across all wrappers of the same atom.  If
  //!          atom's does not already exist, creates it.
  static PyObject * LADA_NAME(getattro)(PyObject *obj, PyObject *name);
  //! \brief sets an attribute.
  //! \details Calls python's _PyObject_GenericGetAttrWithDict(...) with dict
  //!          set to the one contained in atom. This way, we are sure there is
  //!          a single dictionary across all wrappers of the same atom.  If
  //!          atom's does not already exist, creates it.
  static int LADA_NAME(setattro)(PyObject *obj, PyObject *name, PyObject *value);
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
    { Py_INCREF(_self->sequence); return (PyObject*) _self->sequence; }
  static int LADA_NAME(settype)(LADA_TYPE *_self, PyObject *_value, void *closure)
    { to_cpp_sequence_(_value, _self->atom->type); return 0;}
# endif
// Gets dictionary.
static PyObject* LADA_NAME(getdict)(LADA_TYPE *_self, void *_closure)
{
  if(_self->atom->pydict == NULL)
  {
    _self->atom->pydict = PyDict_New();
    if(_self->atom->pydict == NULL) return NULL;
  }
  Py_INCREF(_self->atom->pydict);
  return _self->atom->pydict; 
}
// Sets dictionary.
static int LADA_NAME(setdict)(LADA_TYPE *_self, PyObject *_value, void *_closure)
{
  PyObject* dummy = _self->atom->pydict;
  _self->atom->pydict = _value;
  Py_XINCREF(dummy);
  Py_XDECREF(dummy);
  return 0;
}

static PyObject * LADA_NAME(getattro)(PyObject *obj, PyObject *name)
{
  if(((LADA_TYPE*)obj)->atom->pydict == NULL)
  {
    ((LADA_TYPE*)obj)->atom->pydict = PyDict_New();
    if(((LADA_TYPE*)obj)->atom->pydict == NULL) return NULL;
  }
  return _PyObject_GenericGetAttrWithDict(obj, name, ((LADA_TYPE*)obj)->atom->pydict);
}
static int LADA_NAME(setattro)(PyObject *obj, PyObject *name, PyObject *value)
{
  if(((LADA_TYPE*)obj)->atom->pydict == NULL)
  {
    ((LADA_TYPE*)obj)->atom->pydict = PyDict_New();
    if(((LADA_TYPE*)obj)->atom->pydict == NULL) return -1;
  }
  return _PyObject_GenericSetAttrWithDict(obj, name, value, ((LADA_TYPE*)obj)->atom->pydict);
}

static char LADA_NAME(posdoc)[] = "Position in cartesian coordinates.";
static char LADA_NAME(posname)[] = "pos";
static char LADA_NAME(sitedoc)[] = "Site index (integer).\n\n"
                        "Generally given for a supercell with respect to a backbone lattice."
                        "It is -1 if not initialized, and positive or zero otherwise."
                        "It refers to the original atom in the lattice.";
static char LADA_NAME(dict_doc)[] = "Atom dictionary holding python attributes.";
# if LADA_ATOM_NUMBER  == 0
  static char LADA_NAME(typedoc)[] = "Atomic specie. Must be a string.";
# elif LADA_ATOM_NUMBER == 1
  static char LADA_NAME(typedoc)[] = "List of atomic species.";
# endif 
static char LADA_NAME(sitename)[] = "site";
static char LADA_NAME(freezedoc)[] = "Mask to freeze position or type (unsigned integer).";
static char LADA_NAME(freezename)[] = "freeze";
static char LADA_NAME(type_name)[] = "type";
static char LADA_NAME(dict_name)[] = "__dict__";

extern "C" 
{
  static PyGetSetDef LADA_NAME(getsetters)[] = {
      {LADA_NAME(posname),    (getter)LADA_NAME(getpos), (setter)LADA_NAME(setpos), LADA_NAME(posdoc), NULL},
      {LADA_NAME(sitename),   (getter)LADA_NAME(getsite), (setter)LADA_NAME(setsite), LADA_NAME(sitedoc), NULL},
      {LADA_NAME(freezename), (getter)LADA_NAME(getfreeze),
                              (setter)LADA_NAME(setfreeze), LADA_NAME(freezedoc), NULL}, 
      {LADA_NAME(type_name),  (getter)LADA_NAME(gettype), (setter)LADA_NAME(settype), LADA_NAME(typedoc), NULL}, 
      {LADA_NAME(dict_name),  (getter)LADA_NAME(getdict), (setter)LADA_NAME(setdict), LADA_NAME(typedoc), NULL}, 
      {NULL}  /* Sentinel */
  };
}

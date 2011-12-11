extern "C"
{
  //! Returns cell as a numpy array. 
  static PyObject* LADA_NAME(getcell)(LADA_TYPE *_self, void *closure);
  //! Sets cell from a sequence of 3x3 numbers.
  static int LADA_NAME(setcell)(LADA_TYPE *_self, PyObject *_value, void *_closure);
  // Returns the name of the structure.
  static PyObject* LADA_NAME(getname)(LADA_TYPE *_self, void *closure);
  //! Sets the name of the structure from a string.
  static int LADA_NAME(setname)(LADA_TYPE *_self, PyObject *_value, void *_closure);
  // Returns the energy.
  static PyObject* LADA_NAME(getenergy)(LADA_TYPE *_self, void *closure);
  //! Sets the energy from a number.
  static int LADA_NAME(setenergy)(LADA_TYPE *_self, PyObject *_value, void *_closure);
  // Returns the weight.
  static PyObject* LADA_NAME(getweight)(LADA_TYPE *_self, void *closure);
  //! Sets the weight from a number.
  static int LADA_NAME(setweight)(LADA_TYPE *_self, PyObject *_value, void *_closure);
  // Returns the scale.
  static PyObject* LADA_NAME(getscale)(LADA_TYPE *_self, void *closure);
  //! Sets the scale from a number.
  static int LADA_NAME(setscale)(LADA_TYPE *_self, PyObject *_value, void *_closure);
  // Returns the freeze flag.
  static PyObject* LADA_NAME(getfreeze)(LADA_TYPE *_self, void *closure);
  //! Sets the freeze flag from an unsigned.
  static int LADA_NAME(setfreeze)(LADA_TYPE *_self, PyObject *_value, void *_closure);
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

// Returns cell as a numpy array. 
// Numpy does not implement python's cyclic garbage, hence new wrapper need be
// created each call.
static PyObject* LADA_NAME(getcell)(LADA_TYPE *_self, void *closure)
{
  npy_intp dims[2] = {3, 3};
  int const value = math::numpy::type<math::rMatrix3d::Scalar>::value;
  PyArrayObject* result
    = (PyArrayObject*) PyArray_SimpleNewFromData(2, dims, value, _self->structure->cell.data());
  if(result == NULL) return NULL;
  result->base = (PyObject*)_self;
  Py_INCREF(_self); // Increfed as base of array.
  return (PyObject*)result;
}
// Sets cell from a sequence of three numbers.
static int LADA_NAME(setcell)(LADA_TYPE *_self, PyObject *_value, void *_closure)
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
#   ifdef LADA_NPYITER
#     error LADA_NPYITER is already defined.
#   endif
#   define LADA_NPYITER(TYPE, NUM_TYPE)                                             \
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
          _self->structure->cell(i/3, i%3) = *((TYPE*) PyArray_ITER_DATA(iterator)); \
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
#   undef LADA_NPYITER
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
        while(PyObject* inner = PyObject_GetIter(i_inner))
        {
          if(j >= 3) 
          {
            LADA_PYERROR(TypeError, "Not a 3x3 matrix of numbers.");
            goto inner_error;
          }
          if(PyInt_Check(inner) == 1) _self->structure->cell(i, j) = PyInt_AS_LONG(inner);
          else if(PyFloat_Check(inner) == 1) _self->structure->cell(i, j) = PyInt_AS_LONG(inner);
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
// Returns the name of the structure.
static PyObject* LADA_NAME(getname)(LADA_TYPE *_self, void *closure)
  { return PyString_FromString(_self->structure->name.c_str()); }
// Sets the name of the structure from a string.
static int LADA_NAME(setname)(LADA_TYPE *_self, PyObject *_value, void *_closure)
{
  if(_value == NULL) 
  {
    LADA_PYERROR(TypeError, "Cannot delete name attribute.");
    return -1;
  }
  if(char const * const result = PyString_AsString(_value)) _self->structure->name = result;
  else return -1;
  return 0;
}

// Returns the weight of the structure.
static PyObject* LADA_NAME(getweight)(LADA_TYPE *_self, void *closure)
  { return PyFloat_FromDouble(_self->structure->weight); }
// Sets the weight of the structure from a number.
static int LADA_NAME(setweight)(LADA_TYPE *_self, PyObject *_value, void *_closure)
{
  if(_value == NULL) 
  {
    LADA_PYERROR(TypeError, "Cannot delete name attribute.");
    return -1;
  }
  if(PyFloat_Check(_value)) _self->structure->weight = PyFloat_AS_DOUBLE(_value);
  else if(PyInt_Check(_value)) _self->structure->weight = PyInt_AS_LONG(_value);
  else { LADA_PYERROR(TypeError, "Input is not a number."); return -1; }
  return 0;
}

// Returns the scale of the structure.
static PyObject* LADA_NAME(getscale)(LADA_TYPE *_self, void *closure)
  { return PyFloat_FromDouble(_self->structure->scale); }
// Sets the scale of the structure from a number.
static int LADA_NAME(setscale)(LADA_TYPE *_self, PyObject *_value, void *_closure)
{
  if(_value == NULL) 
  {
    LADA_PYERROR(TypeError, "Cannot delete name attribute.");
    return -1;
  }
  if(PyFloat_Check(_value)) _self->structure->scale = PyFloat_AS_DOUBLE(_value);
  else if(PyInt_Check(_value)) _self->structure->scale = PyInt_AS_LONG(_value);
  else { LADA_PYERROR(TypeError, "Input is not a number."); return -1; }
  return 0;
}
// Returns the energy of the structure.
static PyObject* LADA_NAME(getenergy)(LADA_TYPE *_self, void *closure)
  { return PyFloat_FromDouble(_self->structure->energy); }
// Sets the energy of the structure from a number.
static int LADA_NAME(setenergy)(LADA_TYPE *_self, PyObject *_value, void *_closure)
{
  if(_value == NULL) 
  {
    LADA_PYERROR(TypeError, "Cannot delete name attribute.");
    return -1;
  }
  if(PyFloat_Check(_value)) _self->structure->energy = PyFloat_AS_DOUBLE(_value);
  else if(PyInt_Check(_value)) _self->structure->energy = PyInt_AS_LONG(_value);
  else { LADA_PYERROR(TypeError, "Input is not a number."); return -1; }
  return 0;
}

// Returns the freeze flag.
static PyObject* LADA_NAME(getfreeze)(LADA_TYPE *_self, void *closure)
{
  long result = _self->structure->freeze;
  return PyInt_FromLong(result);
}
// Sets the freeze flag from an unsigned.
static int LADA_NAME(setfreeze)(LADA_TYPE *_self, PyObject *_value, void *_closure)
{
  if(_value == NULL)
  {
    LADA_PYERROR(TypeError, "Cannot delete freeze attribute.");
    return -1;
  }
  long const result = PyInt_AsLong(_value);
  if(result == -1 and PyErr_Occurred() != NULL) return -1;
  if(result < 0)
  {
    PyErr_SetString( PyException<error::ValueError>::exception().ptr(),  
                     "Cannot set freeze to a negative value." );
    return -1;
  }
  _self->structure->freeze = result;
  return 0;
}
// Gets dictionary.
static PyObject* LADA_NAME(getdict)(LADA_TYPE *_self, void *_closure)
{
  if(_self->structure->pydict == NULL)
  {
    _self->structure->pydict = PyDict_New();
    if(_self->structure->pydict == NULL) return NULL;
  }
  Py_INCREF(_self->structure->pydict);
  return _self->structure->pydict; 
}
// Sets dictionary.
static int LADA_NAME(setdict)(LADA_TYPE *_self, PyObject *_value, void *_closure)
{
  Py_XDECREF(_self->structure->pydict);
  _self->structure->pydict = _value;
  Py_XINCREF(_value);
  return 0;
}

static PyObject * LADA_NAME(getattro)(PyObject *obj, PyObject *name)
{
  if(((LADA_TYPE*)obj)->structure->pydict == NULL)
  {
    ((LADA_TYPE*)obj)->structure->pydict = PyDict_New();
    if(((LADA_TYPE*)obj)->structure->pydict == NULL) return NULL;
  }
  return _PyObject_GenericGetAttrWithDict(obj, name, ((LADA_TYPE*)obj)->structure->pydict);
}
static int LADA_NAME(setattro)(PyObject *obj, PyObject *name, PyObject *value)
{
  if(((LADA_TYPE*)obj)->structure->pydict == NULL)
  {
    ((LADA_TYPE*)obj)->structure->pydict = PyDict_New();
    if(((LADA_TYPE*)obj)->structure->pydict == NULL) return -1;
  }
  return _PyObject_GenericSetAttrWithDict(obj, name, value, ((LADA_TYPE*)obj)->structure->pydict);
}

static char LADA_NAME(celldoc)[] = "Position in cartesian coordinates.";
static char LADA_NAME(cellname)[] = "cell";
static char LADA_NAME(weightname)[] = "weight";
static char LADA_NAME(weightdoc)[] = "Weight of the structure, eg for fitting.";
static char LADA_NAME(energyname)[] = "energy";
static char LADA_NAME(energydoc)[] = "energy of the structure, eg for fitting.";
static char LADA_NAME(namename)[] = "name";
static char LADA_NAME(namedoc)[] = "Name of the structure.";
static char LADA_NAME(scalename)[] = "scale";
static char LADA_NAME(scaledoc)[] = "Scale of the structure.\n\n``structure.scale * structure.cell``"
                                    "should be cartesian coordinates in angstrom. Same applies to positions.";
static char LADA_NAME(freezedoc)[] = "Mask to freeze cell coordinate (unsigned integer).";
static char LADA_NAME(freezename)[] = "freeze";
static char LADA_NAME(dictname)[] = "__dict__";
static char LADA_NAME(dictdoc)[] = "Structure dictionary holding python attributes.";

# ifdef LADA_DECLARE 
#   error LADA_DECLARE already defined.
# endif
# define LADA_DECLARE(ATTR)                                          \
    { LADA_NAME(ATTR ## name), (getter)LADA_NAME(get ## ATTR),       \
      (setter)LADA_NAME(set ## ATTR), LADA_NAME(ATTR ## doc), NULL }
extern "C" 
{
  static PyGetSetDef LADA_NAME(getsetters)[] = {
      LADA_DECLARE(cell),
      LADA_DECLARE(freeze),
      LADA_DECLARE(weight),
      LADA_DECLARE(energy),
      LADA_DECLARE(scale),
      LADA_DECLARE(dict), 
      {NULL}
  };
}
# undef LADA_DECLARE

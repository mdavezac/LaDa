extern "C"
{
  //! Returns position as a numpy array. 
  static PyObject* AtomStr_getpos(AtomStr *_self, void *closure);
  //! Sets position from a sequence of three numbers.
  static int AtomStr_setpos(AtomStr *_self, PyObject *_value, void *_closure);
  //! Returns the atomic type. 
  static PyObject* AtomStr_getsite(AtomStr *_self, void *closure);
  //! Sets the atomic type.
  static int AtomStr_setsite(AtomStr *_self, PyObject *_value, void *_closure);
  // Returns the freeze flag.
  static PyObject* AtomStr_getfreeze(AtomStr *_self, void *closure);
  //! Sets the freeze flag from an unsigned.
  static int AtomStr_setfreeze(AtomStr *_self, PyObject *_value, void *_closure);
}

// Returns position as a numpy array. 
static PyObject* AtomStr_getpos(AtomStr *_self, void *closure)
{
  npy_intp dims[1] = {3};
  PyObject *result = PyArray_SimpleNewFromData( 1, dims,
                                                math::numpy::type<math::rVector3d::Scalar>::value,
                                                _self->atom->pos.data());
  if(result != NULL and PyErr_Occurred() == NULL) return result;
  if(result == NULL)
  {
    PyErr_SetString( PyException<error::internal>::exception().ptr(),
                     "Could not create array for position.\n" );
    return NULL;
  }
  Py_DECREF(result);
  return NULL;
}
// Sets position from a sequence of three numbers.
static int AtomStr_setpos(AtomStr *_self, PyObject *_value, void *_closure)
{
  bp::object pos(bp::handle<>(bp::borrowed(_value)));
  if(not is_position(pos)) 
  {
    PyErr_SetString( PyException<error::TypeError>::exception().ptr(),
                     "Input could not be converted to position." );
    return -1;
  }
  extract_position(pos, _self->atom->pos);
  return 0;
}

// Returns the atomic type. 
static PyObject* AtomStr_getsite(AtomStr *_self, void *closure)
  { return PyInt_FromLong(_self->atom->site); }
// Sets the atomic type.
static int AtomStr_setsite(AtomStr *_self, PyObject *_value, void *_closure)
{
  long const result = PyInt_AsLong(_value);
  if(result == -1 and PyErr_Occurred() != NULL) return -1;
  _self->atom->site = result;
  return 0;
}
// Returns the freeze flag.
static PyObject* AtomStr_getfreeze(AtomStr *_self, void *closure)
{
  long result = _self->atom->freeze;
  return PyInt_FromLong(result);
}
// Sets the freeze flag from an unsigned.
static int AtomStr_setfreeze(AtomStr *_self, PyObject *_value, void *_closure)
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

static char posdoc[] = "Position in cartesian coordinates.";
static char posname[] = "pos";
static char sitedoc[] = "Site index (integer).\n\n"
                        "Generally given for a supercell with respect to a backbone lattice."
                        "It is -1 if not initialized, and positive or zero otherwise."
                        "It refers to the original atom in the lattice.";
static char sitename[] = "site";
static char freezedoc[] = "Mask to freeze position or type (unsigned integer).";
static char freezename[] = "freeze";

extern "C" 
{
  static PyGetSetDef AtomStr_getsetters[] = {
      {posname, (getter)AtomStr_getpos, (setter)AtomStr_setpos, posdoc, NULL},
      {sitename, (getter)AtomStr_getsite, (setter)AtomStr_setsite, sitedoc, NULL},
      {freezename, (getter)AtomStr_getfreeze, (setter)AtomStr_setfreeze, freezedoc, NULL}, 
      {NULL}  /* Sentinel */
  };
}

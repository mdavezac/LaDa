extern "C"
{
  //! Returns position as a numpy array. 
  static PyObject* atomstr_getpos(AtomStr *_self, void *closure);
  //! Sets position from a sequence of three numbers.
  static int atomstr_setpos(AtomStr *_self, PyObject *_value, void *_closure);
  //! Returns the atomic type. 
  static PyObject* atomstr_getsite(AtomStr *_self, void *closure);
  //! Sets the atomic type.
  static int atomstr_setsite(AtomStr *_self, PyObject *_value, void *_closure);
  // Returns the freeze flag.
  static PyObject* atomstr_getfreeze(AtomStr *_self, void *closure);
  //! Sets the freeze flag from an unsigned.
  static int atomstr_setfreeze(AtomStr *_self, PyObject *_value, void *_closure);
  // Returns the type of the atom.
  static PyObject* atomstr_gettype(AtomStr *_self, void *closure);
  //! Sets the type from a string.
  static int atomstr_settype(AtomStr *_self, PyObject *_value, void *_closure);
}

// Returns position as a numpy array. 
static PyObject* atomstr_getpos(AtomStr *_self, void *closure)
{
  Py_INCREF(_self->position);
  return (PyObject*)_self->position; 
}
// Sets position from a sequence of three numbers.
static int atomstr_setpos(AtomStr *_self, PyObject *_value, void *_closure)
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
static PyObject* atomstr_getsite(AtomStr *_self, void *closure)
  { return PyInt_FromLong(_self->atom->site); }
// Sets the atomic type.
static int atomstr_setsite(AtomStr *_self, PyObject *_value, void *_closure)
{
  long const result = PyInt_AsLong(_value);
  if(result == -1 and PyErr_Occurred() != NULL) return -1;
  _self->atom->site = result;
  return 0;
}
// Returns the freeze flag.
static PyObject* atomstr_getfreeze(AtomStr *_self, void *closure)
{
  long result = _self->atom->freeze;
  return PyInt_FromLong(result);
}
// Sets the freeze flag from an unsigned.
static int atomstr_setfreeze(AtomStr *_self, PyObject *_value, void *_closure)
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
static PyObject *atomstr_gettype(AtomStr *_self, void *closure)
  { return PyString_FromString(_self->atom->type.c_str()); }
static int atomstr_settype(AtomStr *_self, PyObject *_value, void *closure)
{
  bp::object type(bp::borrowed<>(_value));
  if(not is_specie<std::string>(type))
  {
    PyErr_SetString( PyException<error::TypeError>::exception().ptr(),  
                     "Cannot set type from input. Should be a string." );
    return -1;
  }
  extract_specie<std::string>(type, _self->atom->type);
  return 0;
}

static char posdoc[] = "Position in cartesian coordinates.";
static char posname[] = "pos";
static char sitedoc[] = "Site index (integer).\n\n"
                        "Generally given for a supercell with respect to a backbone lattice."
                        "It is -1 if not initialized, and positive or zero otherwise."
                        "It refers to the original atom in the lattice.";
static char typedoc[] = "Atomic specie as a string.";
static char sitename[] = "site";
static char freezedoc[] = "Mask to freeze position or type (unsigned integer).";
static char freezename[] = "freeze";
static char type_name[] = "type";

extern "C" 
{
  static PyGetSetDef atomstr_getsetters[] = {
      {posname, (getter)atomstr_getpos, (setter)atomstr_setpos, posdoc, NULL},
      {sitename, (getter)atomstr_getsite, (setter)atomstr_setsite, sitedoc, NULL},
      {freezename, (getter)atomstr_getfreeze, (setter)atomstr_setfreeze, freezedoc, NULL}, 
      {type_name, (getter)atomstr_gettype, (setter)atomstr_settype, typedoc, NULL}, 
      {NULL}  /* Sentinel */
  };
}

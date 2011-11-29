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
static PyObject* LADA_NAME(getpos)(LADA_TYPE *_self, void *closure)
{
  Py_INCREF(_self->position);
  return (PyObject*)_self->position; 
}
// Sets position from a sequence of three numbers.
static int LADA_NAME(setpos)(LADA_TYPE *_self, PyObject *_value, void *_closure)
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

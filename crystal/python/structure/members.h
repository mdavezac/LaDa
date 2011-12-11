extern "C" 
{
  //! Returns a representation of the object.
  static PyObject* LADA_NAME(repr)(LADA_TYPE* _self);
  //! If a string atom, returns a sequence, and vice-versa.
  static PyObject* LADA_NAME(cast)(LADA_TYPE* _self, PyObject* _sep);
  //! Returns a deepcopy of the atom.
  static PyObject* LADA_NAME(copy)(LADA_TYPE* _self);
  //! Returns a dictionary with same info as atom.
  static PyObject* LADA_NAME(to_dict)(LADA_TYPE* _self);
  //! Implements deepcopy.
  static PyObject* LADA_NAME(deepcopy)(LADA_TYPE* _self, PyObject* _memo);
  //! Implements shallow copy.
  static PyObject* LADA_NAME(shallowcopy)(LADA_TYPE* _self);
  //! Implements getstate for pickling.
  static PyObject* LADA_NAME(getstate)(LADA_TYPE* _self)
    { return LADA_NAME(to_dict)(_self); }
  //! Implements setstate for pickling.
  static PyObject* LADA_NAME(setstate)(LADA_TYPE* _self, PyObject *_dict);
  //! Implements reduce for pickling.
  static PyObject* LADA_NAME(reduce)(LADA_TYPE* _self);
}

//! Returns a representation of the object.
PyObject* LADA_NAME(repr)(LADA_TYPE* _self)
{
  std::string name(_self->ob_type->tp_name);
  std::ostringstream result;
  // Create structure variable first.
  result << "structure = " << name.substr(name.rfind('.')+1) << "(";
  // Add most attributes in constructor.
  result << _self->structure->scale;
  if(_self->structure->name != "")
    result << ", " << _self->structure->name;
  if(std::abs(_self->structure->weight - 1e0) > 1e-12)
    result << ", weight=" << _self->structure->weight;
  if(std::abs(_self->structure->energy) > 1e-12)
    result << ", energy=" << _self->structure->energy;
  if(_self->structure->freeze != 0)
    result << ", freeze=" << _self->structure->freeze;
  // Including python dynamic attributes.
  if(_self->structure->pydict != NULL)
  {
    if(PyDict_Size(_atom->pydict) > 0)
    {
      PyObject *key, *value;
      Py_ssize_t pos = 0;
      while (PyDict_Next(_atom->pydict, &pos, &key, &value)) 
      {
        PyObject* repr = PyObject_Repr(value);
        if(repr == NULL) return NULL;
        _sstr << ", " << PyString_AsString(key) << "=" << PyString_AsString(repr);
        Py_DECREF(repr);
        if(PyErr_Occurred() != NULL) return NULL;
      }
    }
  }
  // Initialize cell after constructor.
  result << ")\nstructure.cell = [ [" 
         <<         _self->structure->cell(0, 0)
         << ", " << _self->structure->cell(0, 1) 
         << ", " << _self->structure->cell(0, 2) 
         << "],\\\n                   ["
         <<         _self->structure->cell(1, 0)
         << ", " << _self->structure->cell(1, 1) 
         << ", " << _self->structure->cell(1, 2) 
         << "],\\\n                   ["
         <<         _self->structure->cell(2, 0)
         << ", " << _self->structure->cell(2, 1) 
         << ", " << _self->structure->cell(2, 2) 
         << "] ]\n";
# if LADA_ATOM_NUMBER == 0
    typedef LaDa::crystal::StructureData<std::string>::const_iterator t_cit;
# elif LADA_ATOM_NUMBER == 1
    typedef LaDa::crystal::StructureData< std::vector<std::string> >::const_iterator t_cit;
# endif
  // Then add atoms.
  t_cit i_first = _self->structure->begin();
  t_cit const i_end = _self->structure->end();
  result << "structure.add_atom";
  atomstr_repr_impl(*i_first, result);
  for(; i_first != i_end; ++i_first)
  {
    result << "                  ";
#   if LADA_ATOM_NUMBER == 0
      atomstr_repr_impl(*i_first, result);
#   elif LADA_ATOM_NUMBER == 1
      atomsequence_repr_impl(*i_first, result);
#   endif    
    result << "\\\n";
  }
  return PyString_FromString(result.str().c_str());
}

// If a string structure, returns a sequence, and vice-versa.
PyObject *LADA_NAME(cast)(LADA_TYPE *_self, LADA_TYPE *_sep)
{
# if LADA_ATOM_NUMBER == 0
    StructureSequence* result = (StructureSequence*)PyStructureSequence_New();
# elif LADA_ATOM_NUMBER == 1
    StructureStr* result = (StructureStr*)PyStructureStr_New();
# endif
  if(not result) return NULL;
  Py_ssize_t const N(PyTuple_Size(_sep));
  if(N > 1)
  {
    Py_DECREF(result);
    LADA_PYERROR(TypeError, "Input should consist of a single string.");
    return NULL;
  }
  if(N == 1)
  { 
    PyObject *separator = PyTuple_GetItem(_sep, 0);
    char * const sep = PyString_FromString(separator);
    if(sep == NULL) { Py_DECREF(result); return NULL; }
    LaDa::crystal::cast(*_self->structure, result->structure, sep);
  }
  else LaDa::crystal::cast(*_self->structure, result->structure);
  if(PyErr_Occurred() != NULL)
  {
    Py_DECREF(result);
    result = NULL;
  }
  return (PyObject*)result;
}

// Returns a deepcopy of the structure.
PyObject *LADA_NAME(copy)(LADA_TYPE* _self)
{
# if LADA_ATOM_NUMBER == 0
    StructureSequence* result = (StructureSequence*)PyStructureSequence_New();
# elif LADA_ATOM_NUMBER == 1
    StructureStr* result = (StructureStr*)PyStructureStr_New();
# endif
  *result->structure = *_self->structure;
  if(PyErr_Occurred() != NULL) {Py_DECREF(result); return NULL;}
  return (PyObject*)result;
}
// Implements deepcopy of the structure.
PyObject *LADA_NAME(deepcopy)(LADA_TYPE* _self, PyObject* _memo)
{
# if LADA_ATOM_NUMBER == 0
    StructureSequence* result = (StructureSequence*)PyStructureSequence_New();
# elif LADA_ATOM_NUMBER == 1
    StructureStr* result = (StructureStr*)PyStructureStr_New();
# endif
  *result->structure = *_self->structure;
  return (PyObject*)result;
}
// Implements shallow copy.
PyObject *LADA_NAME(shallowcopy)(LADA_TYPE* _self)
{
  Py_INCREF(_self);
  return (PyObject*)_self;
}

// Creates dictionary from atom with shallow copies.
PyObject *LADA_NAME(to_dict)(LADA_TYPE* _self)
{
  PyObject* result = PyDict_New();
  if(not result) return NULL;

# ifdef LADA_ADDITEM 
#   error LADA_ADDITEM already defined.
# endif
# define LADA_ADDITEM(string) \
    item = LADA_NAME(get ## string)((LADA_TYPE*)_self, NULL);           \
    if(item == NULL) goto error;                                        \
    if(PyDict_SetItemString(result, #string, item) < 0) goto erroritem; \
    Py_DECREF(item);
  PyObject* LADA_ADDITEM(cell);
  LADA_ADDITEM(name);
  LADA_ADDITEM(weight);
  LADA_ADDITEM(energy);
  LADA_ADDITEM(scale);
  LADA_ADDITEM(freeze);
#  undef LADA_ADDITEM

  LADA_INNERTYPE::t_Atoms::const_iterator i_atom = _self->structure->begin();
  LADA_INNERTYPE::t_Atoms::const_iterator const i_end = _self->structure->end();
  for(long i(0); i_atom != i_end; ++i_atom, ++i)
  {
    // First, get wrapper to atom.
    // Necessary since items describing the atom will refer to it.
    LADA_ATOM_TYPE* atom = PyAtom_FromAtom(*i_atom);
    if(atom == NULL) { Py_DECREF(result); return NULL; }
    // Then gets dictionary description.
    PyObject *dict = LADA_ATOM_NAME(todict)(atom);
    if(dict == NULL) { Py_DECREF(atom); Py_DECREF(result); return NULL; }
    // Then create pyobject index.
    PyObject *index = PyInt_FromLong(i);
    if(index == NULL) { Py_DECREF(atom); Py_DECREF(result); return NULL; }
    // finally, adds to dictionary.
    bool const error = PyDict_SetItem(result, index, dict) < 0;
    Py_DECREF(atom);
    Py_DECREF(dict);
    Py_DECREF(index);
    if(error) { Py_DECREF(result); return NULL; }
  }
  // Merge attribute dictionary if it exists.
  if(_self->structure->pydict != NULL) PyDict_Merge(result, _self->structure->pydict, 1);

  return result;

  erroritem:
    Py_DECREF(item);
  error: 
    Py_DECREF(result);
    return NULL;
}

// Implements __reduce__ for pickling.
PyObject* LADA_NAME(reduce)(LADA_TYPE* _self)
{
  PyObject *result = PyTuple_New(3);
  if(result == NULL) return NULL;
  Py_INCREF(&LADA_NAME(type));
  if(PyTuple_SET_ITEM(result, 0, (PyObject*)&LADA_NAME(type)) < 0) { Py_DECREF(result); return NULL; }
  PyObject *tuple = PyTuple_New(0);
  if(tuple == NULL) { Py_DECREF(result); return NULL; }
  if(PyTuple_SET_ITEM(result, 1, tuple) < 0) { Py_DECREF(result); return NULL; }
  PyObject *state = LADA_NAME(getstate)(_self);
  if(state == NULL) { Py_DECREF(result); return NULL; }
  if(PyTuple_SET_ITEM(result, 2, state) < 0) { Py_DECREF(result); return NULL; }
  return result;
}

// Implements setstate for pickling.
PyObject* LADA_NAME(setstate)(LADA_TYPE* _self, PyObject *_dict)
{
  PyObject *key, *value;
  Py_ssize_t pos = 0;

  // First find the size of the structure.
  long N(0);
  while (PyDict_Next(_dict, &pos, &key, &value)) 
  {
    if(not PyInt_Check(key)) continue;
    long const i(PyInt_AS_LONG(key));
    if(i > N) N = i;
    else if(i < 0)
    { 
      LADA_PYERROR(internal, "Encountered negative index in dictionary.");
      return NULL;
    }
    if(not PyDict_Check(value)) 
    {
      LADA_PYERROR(internal, "Encounterd integer index which is not an atom dictionary.");
      return NULL;
    }
  }
  // Resize structure.
  try { _self->structure->atoms.resize(N); }
  catch(std::exception &error)
  {
    LADA_PYERROR(internal, ("Encountered error while resizing atoms: " + error.what()).c_str());
    return NULL;
  }
  // some macros to avoid repeating ourselves and to get Str vs Sequence right.
# ifdef LADA_PARSE
#   error LADA_PARSE already exists.
# endif
# define LADA_PARSE(name)                                                            \
    if(index == #name)                                                               \
    {                                                                                \
      if(LADA_NAME(set ## name)(_self, value, NULL) < 0)                             \
      {                                                                              \
        PyErr_Clear();                                                               \
        LADA_PYERROR(ValueError, "Could not unpickle attribute " #name ".");         \
        return NULL;                                                                 \
      }                                                                              \
    }                                                                                
  // Now loop over attributes a second time and sets attributes.
  while (PyDict_Next(_dict, &pos, &key, &value)) 
  {
    if(PyInt_Check(key))
    {
      long const i(PyInt_AS_LONG(key));
      if(i < 0)
      { 
        LADA_PYERROR(internal, "Encountered negative index in dictionary.");
        return NULL;
      }
      LADA_ATOM_TYPE *atom = PyAtom_FromAtom(_self->structure->atoms[i]);
      if(atom == NULL) return NULL;
      if(LADA_ATOM_NAME(setstate)(atom, value) < 0) return NULL;
    }
    else if(PyString_Check(key))
    {
      char const * const index = PyString_AsString(key);
      LADA_PARSE(cell)
      else LADA_PARSE(name)
      else LADA_PARSE(weight)
      else LADA_PARSE(energy)
      else LADA_PARSE(scale)
      else LADA_PARSE(freeze)
      else if(LADA_NAME(setattro)((PyObject*)_self, key, value) < 0) return NULL;
    }
    else if(LADA_NAME(setattro)((PyObject*)_self, key, value) < 0) return NULL;
  }
# undef LADA_PARSE
  Py_RETURN_NONE;
}

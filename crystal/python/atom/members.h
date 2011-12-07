extern "C" 
{
  //! Returns a representation of the object.
  static PyObject* LADA_NAME(repr)(LADA_TYPE* _self);
  //! If a string atom, returns a sequence, and vice-versa.
  static PyObject* LADA_NAME(cast)(LADA_TYPE* _self);
  //! Returns a deepcopy of the atom.
  static PyObject* LADA_NAME(copy)(LADA_TYPE* _self);
  //! Returns a dictionary with same info as atom.
  static PyObject* LADA_NAME(to_dict)(LADA_TYPE* _self);
  //! Implements deepcopy.
  static PyObject* LADA_NAME(deepcopy)(LADA_TYPE* _self, PyObject* _memo);
  //! Implements shallow copy.
  static PyObject* LADA_NAME(shallowcopy)(LADA_TYPE* _self);
}

//! Returns a representation of the object.
PyObject* LADA_NAME(repr)(LADA_TYPE* _self)
{
  std::string name(_self->ob_type->tp_name);
  std::ostringstream result;
  result << name.substr(name.find('.')+1) << "("
         << _self->atom->pos(0) << ", "
         << _self->atom->pos(1) << ", "
         << _self->atom->pos(2);
  if(_self->atom->type.size() > 0) 
#   if LADA_ATOM_NUMBER == 0
      result << ", \'" << _self->atom->type << "\'";
#   elif LADA_ATOM_NUMBER == 1
      {
        std::vector<std::string>::const_iterator i_first = _self->atom->type.begin();
        std::vector<std::string>::const_iterator const i_end = _self->atom->type.end();
        for(; i_first != i_end; ++i_first) result << ", \'" << *i_first << "\'";
      }
#   endif
  if(_self->atom->site != -1)  result << ", site=" << _self->atom->site;
  if(_self->atom->freeze != 0) result << ", freeze=" << _self->atom->freeze;
  if(_self->atom->pydict != NULL)
  {
    if(PyDict_Size(_self->atom->pydict) > 0)
    {
      PyObject *key, *value;
      Py_ssize_t pos = 0;
      while (PyDict_Next(_self->atom->pydict, &pos, &key, &value)) 
      {
        PyObject* repr = PyObject_Repr(value);
        if(repr == NULL) return NULL;
        result << ", " << PyString_AsString(key) << "=" << PyString_AsString(repr);
        Py_DECREF(repr);
        if(PyErr_Occurred() != NULL) return NULL;
      }
    }
  }
  result << ")";
  return PyString_FromString(result.str().c_str());
}


// If a string atom, returns a sequence, and vice-versa.
PyObject *LADA_NAME(cast)(LADA_TYPE *_self)
{
# if LADA_ATOM_NUMBER == 0
    AtomSequence* result = (AtomSequence*)PyAtomSequence_New();
# elif LADA_ATOM_NUMBER == 1
    AtomStr* result = (AtomStr*)PyAtomStr_New();
# endif
  if(not result) return NULL;
  result->atom->pos = _self->atom->pos;
  result->atom->freeze = _self->atom->freeze;
  result->atom->site = _self->atom->site;
  if(_self->atom->type.size() > 0)
#   if LADA_ATOM_NUMBER == 1
      result->atom->type = LaDa::crystal::details::print_occupation(_self->atom->type);
#   elif LADA_ATOM_NUMBER == 0
      result->atom->type.push_back(_self->atom->type);
#   endif
  if(_self->atom->pydict != NULL)
  {
    if(PyDict_Size(_self->atom->pydict) > 0)
    {  
      Py_XDECREF(result->atom->pydict);
      result->atom->pydict = NULL;
      PyObject* copymod = PyImport_ImportModule("copy");
      if(copymod == NULL) { Py_DECREF(result); return NULL; }
      PyObject *deepcopystr = PyString_FromString("deepcopy");
      if(not deepcopystr) { Py_DECREF(result); Py_DECREF(copymod); return NULL; }
      result->atom->pydict = PyObject_CallMethodObjArgs(copymod, deepcopystr, _self->atom->pydict, NULL);
      Py_DECREF(copymod);
      Py_DECREF(deepcopystr);
      if(not result->atom->pydict) { Py_DECREF(result); return NULL; }
    }
  }
  return (PyObject*)result;
}

// Returns a deepcopy of the atom.
PyObject *LADA_NAME(copy)(LADA_TYPE* _self)
{
# if LADA_ATOM_NUMBER == 0
    LADA_TYPE* result = (LADA_TYPE*)PyAtomStr_New();
# elif LADA_ATOM_NUMBER == 1
    LADA_TYPE* result = (LADA_TYPE*)PyAtomSequence_New();
# endif
  *result->atom = *_self->atom;
  if(PyErr_Occurred() != NULL) {Py_DECREF(result); return NULL;}
  return (PyObject*)result;
}
// Implements deepcopy of the atom.
PyObject *LADA_NAME(deepcopy)(LADA_TYPE* _self, PyObject* _memo)
{
# if LADA_ATOM_NUMBER == 0
    LADA_TYPE* result = (LADA_TYPE*)PyAtomStr_New();
# elif LADA_ATOM_NUMBER == 1
    LADA_TYPE* result = (LADA_TYPE*)PyAtomSequence_New();
# endif
  result->atom->pos = _self->atom->pos;
  result->atom->type = _self->atom->type;
  result->atom->freeze = _self->atom->freeze;
  result->atom->site = _self->atom->site;
  if(_self->atom->pydict != NULL)
  {
    Py_XDECREF(result->atom->pydict);
    result->atom->pydict = NULL;
    PyObject* copymod = PyImport_ImportModule("copy");
    if(copymod == NULL) { Py_DECREF(result); return NULL; }
    PyObject *deepcopystr = PyString_FromString("deepcopy");
    if(deepcopystr == NULL) { Py_DECREF(result); Py_DECREF(copymod); return NULL; }
    result->atom->pydict = PyObject_CallMethodObjArgs(copymod, deepcopystr, _self->atom->pydict, _memo, NULL);
    Py_DECREF(copymod);
    Py_DECREF(deepcopystr);
    if(result->atom->pydict == NULL) { Py_DECREF(result); return NULL; }
  }
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
  PyObject* LADA_ADDITEM(pos);
  LADA_ADDITEM(type);
  LADA_ADDITEM(site);
  LADA_ADDITEM(freeze);
# undef LADA_ADDITEM

  // Merge attribute dictionary if it exists.
  if(_self->atom->pydict != NULL) PyDict_Merge(result, _self->atom->pydict, 1);

  return result;

  erroritem:
    Py_DECREF(item);
  error: 
    Py_DECREF(result);
    return NULL;
}

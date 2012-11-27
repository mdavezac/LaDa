//! Implements atom representation. 
int lada_atom_repr_impl(PyAtomObject const &_atom, std::ostringstream &_sstr);
//! Returns a representation of the object.
PyObject* lada_atom_repr(PyAtomObject* _self)
{
  std::string name(_self->ob_type->tp_name);
  std::ostringstream sstr;
  sstr << name.substr(name.rfind('.')+1);
  if(lada_atom_repr_impl(*_self, sstr) < 0) return NULL;
  return PyString_FromString(sstr.str().c_str());
}
//! Returns a deepcopy of the atom.
PyObject* lada_atom_copy(PyAtomObject* _self)
  { return (PyObject*) copy_atom(_self, NULL); }
//! Implements shallow copy.
PyObject* lada_atom_shallowcopy(PyAtomObject* _self)
  { Py_INCREF(_self); return (PyObject*)_self; }
//! Returns a dictionary with same info as atom.
PyObject* lada_atom_to_dict(PyAtomObject* _self);
//! Implements getstate for pickling.
PyObject* lada_atom_getstate(PyAtomObject* _self);
//! Implements setstate for pickling.
PyObject* lada_atom_setstate(PyAtomObject* _self, PyObject *_dict);
//! Implements reduce for pickling.
PyObject* lada_atom_reduce(PyAtomObject* _self);


// Implements atom representation. 
int lada_atom_repr_impl(PyAtomObject const &_atom, std::ostringstream &_sstr)
{
  _sstr << "("
         << _atom.pos(0) << ", "
         << _atom.pos(1) << ", "
         << _atom.pos(2);
  if(_atom.type != Py_None)
  {
    PyObject* repr = PyObject_Repr(_atom.type);
    if(repr == NULL) return -1;
    _sstr << ", " << PyString_AsString(repr);
  }
  if(_atom.pydict != NULL)
  {
    if(PyDict_Size(_atom.pydict) > 0)
    {
      PyObject *key, *value;
      Py_ssize_t pos = 0;
      while (PyDict_Next(_atom.pydict, &pos, &key, &value)) 
      {
        PyObject* repr = PyObject_Repr(value);
        if(repr == NULL) return -1;
        _sstr << ", " << PyString_AsString(key) << "=" << PyString_AsString(repr);
        Py_DECREF(repr);
        if(PyErr_Occurred() != NULL) return -1;
      }
    }
  }
  _sstr << ")";
  return 0;
}

// Creates dictionary from atom with shallow copies.
PyObject *lada_atom_to_dict(PyAtomObject* _self)
{
  PyObject* result = PyDict_New();
  if(not result) return NULL;

  PyObject *item = lada_atom_getpos((PyAtomObject*)_self, NULL);
  if(item == NULL) { Py_DECREF(result); return NULL; }
  bool const iserror = PyDict_SetItemString(result, "pos", item) < 0; 
  Py_DECREF(item);
  if(iserror) goto error;

  // Merge attribute dictionary if it exists.
  if(_self->pydict != NULL)
    if(PyDict_Merge(result, _self->pydict, 1) < 0) goto error;
  if(PyDict_SetItemString(result, "type", _self->type) < 0) goto error;

  return result;

  error: 
    Py_DECREF(result);
    return NULL;
}

// Implements __reduce__ for pickling.
PyObject* lada_atom_reduce(PyAtomObject* _self)
{
  // Creates return tuple of three elements.
  PyObject * const result = PyTuple_New(3);
  if(result == NULL) return NULL;
  // first element is the object type.
  if(PyTuple_SET_ITEM(result, 0, PyObject_Type((PyObject*)_self)) < 0)
    { Py_DECREF(result); return NULL; }
  // Second element is a null tuple, argument to the callable type above.
  PyObject * const tuple = PyTuple_New(0);
  if(tuple == NULL) { Py_DECREF(result); return NULL; }
  if(PyTuple_SET_ITEM(result, 1, tuple) < 0) { Py_DECREF(result); return NULL; }
  // Third element is the state of this object.
  char getstate[] = "__getstate__";
  PyObject * state = PyObject_CallMethod((PyObject*)_self, getstate, NULL, NULL);
  if(state == NULL) { Py_DECREF(result); return NULL; }
  if(PyTuple_SET_ITEM(result, 2, state) < 0) { Py_DECREF(result); return NULL; }
  // Finally, returns.
  return result;
}
// Implements getstate for pickling.
PyObject* lada_atom_getstate(PyAtomObject* _self)
{
  PyObject * const pos = lada_atom_getpos(_self, NULL);
  if(pos == NULL) return NULL;
  PyObject * const type = lada_atom_gettype(_self, NULL);
  if(type == NULL) { Py_DECREF(pos); return NULL; }
  PyObject * dict = _self->pydict == NULL ? Py_None: _self->pydict; 
  Py_INCREF(dict);

  PyObject *const result = PyTuple_New(3);
  if(result == NULL) goto error0;
  if(PyTuple_SET_ITEM(result, 0, pos) < 0) goto error1;
  if(PyTuple_SET_ITEM(result, 1, type) < 0) goto error2;
  if(PyTuple_SET_ITEM(result, 2, dict) < 0) goto error3;
  return result;

  error0:
    Py_DECREF(pos);
  error1:
    Py_DECREF(type);
  error2:
    Py_DECREF(dict);
  error3:
    Py_DECREF(result);
    return NULL;
}

// Implements setstate for pickling.
PyObject* lada_atom_setstate(PyAtomObject* _self, PyObject *_tuple)
{
  if(not PyTuple_Check(_tuple))
  {
    LADA_PYERROR(TypeError, "Atom.__setstate__: expected a tuple.");
    return NULL;
  }
  if(PyTuple_Size(_tuple) != 3)
  {
    LADA_PYERROR(TypeError, "Atom.__setstate__: expected a tuple of three elements.");
    return NULL;
  }
  if(PyObject * const item = PyTuple_GET_ITEM(_tuple, 0))
    { lada_atom_setpos(_self, item, 0); }
  else return NULL;
  if(PyObject * const item = PyTuple_GET_ITEM(_tuple, 1))
    { lada_atom_settype(_self, item, 0); }
  else return NULL;
  if(PyObject * const item = PyTuple_GET_ITEM(_tuple, 2))
  {
    if(item != Py_None)
    {
      PyObject *dummy = _self->pydict;
      _self->pydict = item;
      Py_INCREF(_self->pydict);
      Py_DECREF(dummy);
    }
  }
  else return NULL;
  Py_RETURN_NONE;
}

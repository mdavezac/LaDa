namespace LaDa
{
  namespace crystal
  {
    extern "C" 
    { 
      //! Function to deallocate a string atom.
      static void lada_atom_dealloc(PyAtomObject *_self);
      //! Function to initialize a string atom.
      static int lada_atom_init(PyAtomObject* _self, PyObject* _args, PyObject *_kwargs);
      //! Traverses to back-reference.
      static int lada_atom_traverse(PyAtomObject *_self, visitproc visit, void *arg)
        { Py_VISIT(_self->pydict); Py_VISIT(_self->type); return 0; }
      //! Clears back reference.
      static int lada_atom_gcclear(PyAtomObject *self) { Py_CLEAR(self->pydict); Py_CLEAR(self->type); return 0; }
    }
 
    // Function to deallocate a string atom.
    static void lada_atom_dealloc(PyAtomObject *_self)
    {
      if(_self->weakreflist != NULL)
        PyObject_ClearWeakRefs((PyObject *) _self);
     
      lada_atom_gcclear(_self);

      // Calls c++ destructor explicitely.
      PyTypeObject *ob_type = _self->ob_type;
      _self->~PyAtomObject();

      ob_type->tp_free((PyObject*)_self);
    }
 
    // Function to initialize a string atom.
    static int lada_atom_init(PyAtomObject* _self, PyObject* _args, PyObject *_kwargs)
    {
      int found_position = 0;
      bool found_type = false;
      Py_ssize_t const N = PyTuple_Size(_args);
 
      if(N > 0) 
      {
        if(N < 3)
        {
          LADA_PYERROR(TypeError, "Atom(...): Expect at least three arguments.");
          return -1;
        }
        for(Py_ssize_t i(0); i < 3; ++i)
        {
          PyObject *item = PyTuple_GET_ITEM(_args, i);
          if(PyInt_Check(item) == 1)  _self->pos[i] = PyInt_AS_LONG(item);
          else if(PyFloat_Check(item) == 1)  _self->pos[i] = PyFloat_AS_DOUBLE(item);
          else
          {
            LADA_PYERROR( TypeError,
                          "Atom(...): Expects the first three arguments to be the position. "
                          "Or, give everything as keywords" );
            return -1;
          }
        }
        if(N == 4) 
        {
          // deletes current type.
          PyObject* dummy = _self->type;
          _self->type = PyTuple_GET_ITEM(_args, 3);
          Py_INCREF(_self->type);
          Py_DECREF(dummy);
        }
        else if(N > 4)
        {
          // deletes current type.
          PyObject* dummy = _self->type;
          PyObject *slice = PyTuple_GetSlice(_args, 3, N);
          if(slice == NULL) return -1;
          _self->type = PyList_New(N-3);
          if(_self->type == NULL)
          {
            _self->type = dummy;
            return -1;
          }
          Py_DECREF(dummy);
          if(PyList_SetSlice(_self->type, 0, N-3, slice) < -1) return -1;
        }
      }
 
      // check for keywords.
      if(_kwargs == NULL) return 0;
 
      // Sanity check first.
      if( N >= 3 and PyDict_GetItemString(_kwargs, "pos") != NULL)
      {
        LADA_PYERROR(TypeError, "Atom(...): Cannot set position both from keyword and argument.");
        return -1;
      }
      // Sanity check first.
      if( N >= 4 and PyDict_GetItemString(_kwargs, "type") != NULL)
      {
        LADA_PYERROR(TypeError, "Atom(...): Cannot set type both from keyword and argument.");
        return -1;
      }
 
      PyObject *key, *value;
      Py_ssize_t pos = 0;
      while (PyDict_Next(_kwargs, &pos, &key, &value)) 
        if(PyObject_SetAttr((PyObject*)_self, key, value) < 0) return -1;
      return 0;
    }
  } // namespace Crystal
} // namespace LaDa

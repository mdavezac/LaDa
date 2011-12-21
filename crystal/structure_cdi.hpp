namespace LaDa
{
  namespace crystal
  {
    extern "C" 
    { 
      //! Function to deallocate a string atom.
      static void structure_dealloc(StructureData *_self);
      //! Function to initialize a string atom.
      static int structure_init(StructureData* _self, PyObject* _args, PyObject *_kwargs);
      //! Traverses to back-reference.
      static int structure_traverse(StructureData *_self, visitproc _visit, void *_arg);
      //! Clears back reference.
      static int structure_gcclear(StructureData *_self);
    }
  
    // Function to deallocate a string atom.
    static void structure_dealloc(StructureData *_self)
    {
      if(_self->weakreflist != NULL)
        PyObject_ClearWeakRefs((PyObject *) _self);
     
      structure_gcclear(_self);
  
      _self->ob_type->tp_free((PyObject*)_self);
    }
  
    // Function to initialize an atom.
    static int structure_init(StructureData* _self, PyObject* _args, PyObject *_kwargs)
    {
      Py_ssize_t const N = PyTuple_Size(_args);
      if(N > 1 and N != 9)
      {
        LADA_PYERROR(TypeError, "Unexpected argument: arguments should represent the cell "
                                "and be given as a matrix, or as a series of 9 numbers." );
      }
      else if(N == 1)
      {
        if(_kwargs != NULL and PyDict_GetItemString(_kwargs, "cell") != NULL)
        {
          LADA_PYERROR(TypeError, "Cell given as both argument and keyword.");
          return -1;
        }
        PyObject *item = PyTuple_GetItem(_args, 0);
        if(item == NULL) return -1;
        if(structure_setcell(_self, item, NULL) < 0) return -1;
      }
      else if(N == 9 and structure_setcell(_self, _args, NULL) < 0) return -1;

      if(_kwargs == NULL) return 0;
      PyObject *key, *value;
      Py_ssize_t pos = 0;
      while (PyDict_Next(_kwargs, &pos, &key, &value)) 
        if(PyObject_SetAttr((PyObject*)_self, key, value) < 0) return -1;
      return 0;
    }
  
    static int structure_traverse(StructureData *self, visitproc visit, void *arg)
    {
      Py_VISIT(self->pydict);
      return 0;
    }
  
    static int structure_gcclear(StructureData *self)
    { 
      Py_CLEAR(self->pydict);
      return 0;
    }
  }
}

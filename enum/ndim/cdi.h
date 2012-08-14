namespace LaDa
{
  namespace enumeration
  {
    extern "C" 
    { 
      //! Function to deallocate a string atom.
      static void ndimiterator_dealloc(NDimIterator *_self);
      //! Function to initialize a string atom.
      static int ndimiterator_init(NDimIterator* _self, PyObject* _args, PyObject *_kwargs);
      //! Traverses to back-reference.
      static int ndimiterator_traverse(NDimIterator *_self, visitproc _visit, void *_arg);
      //! Clears back reference.
      static int ndimiterator_gcclear(NDimIterator *_self);
    }
  
    // Function to deallocate a string atom.
    static void ndimiterator_dealloc(NDimIterator *_self)
    {
      ndimiterator_gcclear(_self);
  
      // calls destructor explicitely.
      PyTypeObject* ob_type = _self->ob_type;
      _self->~NDimIterator();

      ob_type->tp_free((PyObject*)_self);
    }
  
    // Function to initialize an atom.
    static int ndimiterator_init(NDimIterator* _self, PyObject* _args, PyObject *_kwargs)
    {
      if(_args == NULL)
      {
        LADA_PYERROR(TypeError, "NDimIterator expects at least one argument.");
        return -1;
      }
      Py_ssize_t const N = PyTuple_Size(_args);
      if(N == 0)
      {
        LADA_PYERROR(TypeError, "NDimIterator expects at least one argument.");
        return -1;
      }
      if(_kwargs != NULL and PyDict_Size(_kwargs))
      {
        LADA_PYERROR(TypeError, "NDimIterator does not expect keyword arguments.");
        return -1;
      }
      for(Py_ssize_t i(0); i < N; ++i)
      {
        PyObject *item = PyTuple_GET_ITEM(_args, i);
        if(PyLong_Check(item)) _self->ends.push_back( PyLong_AS_LONG(item));
        else if(PyInt_Check(item)) _self->ends.push_back( PyInt_AS_LONG(item));
        else
        {
          LADA_PYERROR(TypeError, "Unknown type in NDimIterator.");
          return -1;
        }
        if(_self.ends.back() < 0)
        {
          LADA_PYERROR(TypeError, "NDimIterator does not expect negative arguments.");
          return -1;
        }
        if(_self.ends.back() > 5)
        {
          LADA_PYERROR(TypeError, "Argument larger than 5. "
                       "Recompile LaDa if that is what you want.");
          return -1;
        }
      } 
      _self->vector.reserve(5);
      std::fill(_self->vector.begin(), _self->vector.end(), 1);
      return 0;
    }
  
    static int ndimiterator_traverse(NDimIterator *self, visitproc visit, void *arg)
    {
      Py_VISIT(self->yielded);
      return 0;
    }
  
    static int ndimiterator_gcclear(NDimIterator *self)
    { 
      Py_CLEAR(self->yielded);
      return 0;
    }
  }
}

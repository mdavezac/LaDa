extern "C" 
{ 
  //! Function to deallocate a string atom.
  static void LADA_NAME(dealloc)(LADA_TYPE *_self);
  //! Function to allocate a string atom.
  static PyObject* LADA_NAME(new)(PyTypeObject *_type, PyObject *_args, PyObject *_kwargs);
  //! Function to initialize a string atom.
  static int LADA_NAME(init)(LADA_TYPE* _self, PyObject* _args, PyObject *_kwargs);
  //! Traverses to back-reference.
  static int LADA_NAME(traverse)(LADA_TYPE *_self, visitproc _visit, void *_arg);
  //! Clears back reference.
  static int LADA_NAME(gcclear)(LADA_TYPE *_self);
}

// Function to deallocate a string atom.
static void LADA_NAME(dealloc)(LADA_TYPE *_self)
{
  if(_self->weakreflist != NULL)
    PyObject_ClearWeakRefs((PyObject *) _self);

  LADA_NAME(gcclear)(_self);
  boost::shared_ptr< LaDa::crystal::AtomData< std::string > > dummy;
  dummy.swap(_self->atom);
  _self->ob_type->tp_free((PyObject*)_self);
}

//! Function to allocate a string atom.
static PyObject* LADA_NAME(new)(PyTypeObject *_type, PyObject *_args, PyObject *_kwargs)
{
  PyObject *const pydict = PyDict_New();
  if(pydict == NULL) return NULL;
  LADA_TYPE *self;
  self = (LADA_TYPE*)_type->tp_alloc(_type, 0);
  if(self == NULL) return NULL;
  self->weakreflist = NULL;
  boost::shared_ptr< crystal::AtomData<std::string> > dummy(new LaDa::crystal::AtomData<std::string>);
  if(not dummy)
  {
    Py_DECREF(self);
    LADA_PYERROR(internal, "Could not create atom.\n" );
    return NULL;
  } 
  self->atom.swap(dummy);
  self->dictionary = pydict;


# if LADA_ATOM_NUMBER == 1
    self->sequence = (Sequence*)sequence_type.tp_alloc(sequence_type, 0);
    if(self->sequence == NULL)
    {
      Py_DECREF(self->position);
      Py_DECREF(self);
      return NULL;
    }
    self->sequence->ptr_seq = &self->atom->type;
    self->sequence->ptr_base = self;
    Py_INCREF(self); // increfed because of sequence base.
# endif

  return (PyObject*) self;
}

// Function to initialize a string atom.
static int LADA_NAME(init)(LADA_TYPE* _self, PyObject* _args, PyObject *_kwargs)
{
  int found_position = 0;
  bool found_type = false;
  Py_ssize_t const N = PyTuple_Size(_args);
  // loop over arguments. 
  // Try and identify whether arguments are positions, or occupation string.
  // Leaves a certain amount of lattitude, in that both positions can be given
  // a sequences, or each component given one at a time. Makes for more complex
  // unpacking though.
  for(Py_ssize_t i(0); i < N; ++i)
  {
    PyObject *const item = PyTuple_GET_ITEM(_args, i);
    if(PyString_Check(item) == 1)
    {
#     if LADA_ATOM_NUMBER == 0
        if(found_type) { LADA_PYERROR(TypeError, "More than one occupation given."); return -1; }
        _self->atom->type = PyString_AS_STRING(item);
#     elif LADA_ATOM_NUMBER == 1
        _self->atom->type.push_back(PyString_AS_STRING(item));
#     endif
      found_type = true;
      continue;
    }
    else if(i < 3 and found_position < 3 and found_type == false)
    {
      if(PyInt_Check(item) == 1) 
      {
        if(found_position != i)
        {
          LADA_PYERROR(TypeError, "Position can be given either as sequence or as three arguments.");
          return -1;
        }
        _self->atom->pos[found_position] = PyInt_AS_LONG(item);
        ++found_position;
        continue;
      }
      else if(PyFloat_Check(item) == 1)
      {
        if(found_position != i)
        {
          LADA_PYERROR(TypeError, "Position can be given either as sequence or as three arguments.");
          return -1;
        }
        _self->atom->pos[found_position] = PyFloat_AS_DOUBLE(item);
        ++found_position;
        continue;
      }
    }
    if(PyObject *iterator = PyObject_GetIter(item))
    {
      size_t j(0);
      bool is_pos = false;
      bool is_type = false;
      while(PyObject* inner_item = PyIter_Next(iterator))
      {
        if(i == 0 and PyFloat_Check(inner_item) == 1 or PyInt_Check(inner_item) == 1)
        {
          if(is_type == true) 
          {
            LADA_PYERROR(TypeError, "Cannot tell wether tuple/sequence argument is type or position.");
            goto error;
          }
          if(j != found_position)
          {
            LADA_PYERROR(TypeError, "Tuple/sequence argument defining position should "
                                    "contain exactly three separate argument.");
            goto error;
          }
          if(found_position >= 3)
          {
            LADA_PYERROR(TypeError, "Tuple/sequence argument defining position should "
                                    "contain exactly three separate argument.");
            goto error;
          }
          is_pos == true;
          if(PyInt_Check(inner_item) == 1) _self->atom->pos[j] = PyInt_AS_LONG(inner_item); 
          else if(PyFloat_Check(inner_item) == 1) _self->atom->pos[j] = PyFloat_AS_DOUBLE(inner_item); 
          else 
          {
            LADA_PYERROR(TypeError, "Found sequence specifying position, "
                                    "but could not make sense of its components as integers or floats.");
            goto error;
          }
          ++found_position;
        } // if integer or float and correct position for pos.
#       if LADA_ATOM_NUMBER == 1         
          else if(PyString_Check(inner_item) == 1)
          {
            if(j == 0 and found_type)
            {
              LADA_PYERROR(TypeError, "Unexpected extra arguments.");
              goto error;
            }
            if(char * const string = PyString_AsString(inner_item))
            {
              _self->atom->type.push_back(inner_item);
              found_type = true;
            }
            else
            {
              LADA_PYERROR(TypeError, "Found sequence specifying position, "
                                      "but could not make sense of its components as integers or floats.");
              goto error;
            }
          } // if list of types.
#       endif
        else 
        {
          LADA_PYERROR(TypeError, "Found a sequence, but not the sense of it.");
          return -1;
        } // don't know what to make of input.
        ++j;
        Py_DECREF(inner_item);
        continue;
        error:
          Py_DECREF(inner_item);
          Py_DECREF(iterator);
          return -1;
      } // loop over inner tuple.
      Py_DECREF(iterator);
    } // end of iterable sequence.
    else 
    {
      LADA_PYERROR(TypeError, "Could not make sense of arguments.");
      return -1;
    }
  } // loop over arguments.
  if(found_position < 3 and found_position > 0)
  {
    LADA_PYERROR(TypeError, "Could not find all positions in arguments.");
    return -1;
  }
  if(_kwargs == NULL) return 1;
  
  if(PyObject *type = PyDict_GetItemString(_kwargs, "type"))
  {
    if(found_type)
    {
      LADA_PYERROR(TypeError, "Type given both as argument and keyword argument.");
      return -1;
    }
#   if LADA_ATOM_NUMBER == 0
      if(PyString_Check(type) == 1) _self->atom->type = PyString_AS_STRING(type);
#   elif LADA_ATOM_NUMBER == 1
      if(PyObject* return_ = to_cpp_sequence_(input, self->atom->type)) Py_DECREF(return_);
#   endif
    else 
    {
      LADA_PYERROR(TypeError, "Could not make sense of type keyword argument.");
      return -1;
    }
    PyDict_DelItemString(_kwargs, "type");
  }
  if(PyObject *pos = PyDict_GetItemString(_kwargs, "pos"))
  {
    if(found_position != 0)
    {
      LADA_PYERROR(TypeError, "Found position as argument and keyword.");
      return -1;
    }
    if(PyObject *iterator = PyObject_GetIter(pos))
    {
      size_t j(0);
      while(PyObject *item = PyIter_Next(iterator))
      {
        if(j >= 3)
        {
          LADA_PYERROR(TypeError, "More than three components given for position.");
          Py_DECREF(item);
          Py_DECREF(iterator);
          return -1;
        }
        if(PyInt_Check(item) == 1) _self->atom->pos[j] = PyInt_AS_LONG(item); 
        else if(PyFloat_Check(item) == 1) _self->atom->pos[j] = PyFloat_AS_DOUBLE(item);
        else 
        {
          Py_DECREF(item); Py_DECREF(iterator); 
          LADA_PYERROR(TypeError, "pos keyword argument should be a sequence of three numbers.");
          return -1;
        }
        Py_DECREF(item);
        ++j;
      }
      Py_DECREF(iterator);
      found_position = 3;
    }
    else
    {
      LADA_PYERROR(TypeError, "pos keyword argument should be a sequence of three numbers.");
      return -1;
    }
    PyDict_DelItemString(_kwargs, "pos");
  }
  if(PyObject *site = PyDict_GetItemString(_kwargs, "site"))
  {
    if(PyInt_Check(site) == 0)
    {
      LADA_PYERROR(TypeError, "site keyword argument should be an integer.");
      return -1;
    }
    _self->atom->site = PyInt_AS_LONG(site);
    PyDict_DelItemString(_kwargs, "site");
  }
  if(PyObject *freeze = PyDict_GetItemString(_kwargs, "freeze"))
  {
    if(PyInt_Check(freeze) == 0)
    {
      LADA_PYERROR(TypeError, "freeze keyword argument should be an integer.");
      return -1;
    }
    _self->atom->freeze = PyInt_AS_LONG(freeze);
    PyDict_DelItemString(_kwargs, "freeze");
  }
  // Now additional attributes.
  int const result = PyDict_Merge(_self->dictionary, _kwargs, 1);
  return result;
}

static int LADA_NAME(traverse)(LADA_TYPE *self, visitproc visit, void *arg)
{
  Py_VISIT(self->dictionary);
# if LADA_ATOM_NUMBER == 1
    Py_VISIT(self->sequence);
# endif
  return 0;
}

static int LADA_NAME(gcclear)(LADA_TYPE *self)
{ 
  Py_CLEAR(self->dictionary);
# if LADA_ATOM_NUMBER == 1
    Py_CLEAR(self->sequence);
# endif
  return 0;
}

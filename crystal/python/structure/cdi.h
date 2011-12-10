extern "C" 
{ 
  //! Function to deallocate a string atom.
  static void LADA_NAME(dealloc)(LADA_TYPE *_self);
  //! Function to allocate a string atom.
  static PyObject* LADA_NAME(new)(PyTypeObject *_type, PyObject *_args, PyObject *_kwargs)
#    if LADA_ATOM_NUMBER == 0
       { return PyStructureStr_New(); }
#    elif LADA_ATOM_NUMBER == 1
       { return PyStructureSequence_New(); }
#    endif
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
 
  // Reference to wrapper in the wrapped object is not owned by the wrapped
  // object. It is set to NULL when the wrapper is destroyed. 
  if(_self->structure) _self->structure->pyself = NULL;

  LADA_NAME(gcclear)(_self);

  _self->ob_type->tp_free((PyObject*)_self);
}

// Creates a new atom and its wrapper.
#if LADA_ATOM_NUMBER == 0
  PyObject* PyStructureStr_New()
#elif LADA_ATOM_NUMBER == 1
  PyObject* PyStructureSequence_New()
#endif
{
  LADA_TYPE* result = (LADA_TYPE*)LADA_NAME(type).tp_alloc(&LADA_NAME(type), 0);
  if(result == NULL) return NULL;
  
  // set everything to null, just in case we exit to fast.
  result->weakreflist = NULL;
  // Now starts setting things up.
# if LADA_ATOM_NUMBER == 0
    typedef LaDa::crystal::StructureData<std::string> t_Structure;
# elif LADA_ATOM_NUMBER == 1
    typedef LaDa::crystal::StructureData< std::vector<std::string> > t_Structure;
# endif
  result->structure.reset(new(std::nothrow) t_Structure);
  if(not result->structure)
  {
    Py_DECREF(result);
    LADA_PYERROR(internal, "Could not create atom.\n" );
    return NULL;
  }
  // Reference to wrapper in the wrapped object is not owned by the wrapped
  // object. It is set to NULL when the wrapper is destroyed. 
  result->structure->pyself = (PyObject*)result;
  return (PyObject*) result;
}

// Function to initialize a string atom.
static int LADA_NAME(init)(LADA_TYPE* _self, PyObject* _args, PyObject *_kwargs)
{
  Py_ssize_t const N = PyTuple_Size(_args);
  Py_ssize_t i(0);

# ifdef LADA_SET
#  error LADA_SET already defined
# endif
# if LADA_SET(attr, default_)                                              \
    {                                                                      \
      /* Check for keyword argument first. */                              \
      PyObject* item = PyDict_GetItemString(_kwargs, #attr);               \
      bool const iskw = item != NULL;                                      \
      /* If not a keyword, checked for next unparsed argument. */          \
      if(N > i and not iskw) cell = PyTuple_GetItem(_args, i);             \
      if(item != NULL)                                                     \
      {                                                                    \
        /* Check if we can set the object using standard get/set method */ \
        if(LADA_NAME(set ## attr)(_self, item, NULL) == -1)                \
        {                                                                  \
          /* If a keyword and could not parse, then return error. */       \
          if(iskw) return -1;                                              \
          /* If this came from an argument, then could be the argument     \
           * wasn't meant for this attribute but the next.                 \
           * Clear error and set default. */                               \
          PyErr_Clear();                                                   \
          _self->structure->attr = default_;                               \
        }                                                                  \
        /* Increment unparsed argument position */                         \
        if(not iskw) ++i;                                                  \
        /* Remove keyword argument to make adding dynamic python attribute \
         * easier later on. */                                             \
        else if(PyDict_DelItemString(_kwargs, #attr) == -1) return -1;     \
      }                                                                    \
      /* Neither keyword nor argument: set to default. */                  \
      else _self->structure->attr = default_;                              \
    }
  LADA_SET(cell,   LaDa::math::rMatrix3d::Identity());
  LADA_SET(scale,  1e0);
  LADA_SET(name,    "");
  LADA_SET(weight, 1e0);
  LADA_SET(energy, 0e0);
  LADA_SET(freeze,   0);
# undef LADA_SET
  // Now additional attributes.
  if(_kwargs == NULL) return 0;
  _self->structure->pydict = PyDict_New();
  return PyDict_Merge(_self->atom->pydict, _kwargs, 1);
}

static int LADA_NAME(traverse)(LADA_TYPE *self, visitproc visit, void *arg)
{
  Py_VISIT(self->structure->pydict);
  return 0;
}

static int LADA_NAME(gcclear)(LADA_TYPE *self)
{ 
  // Following line basically calls Py_CLEAR if atom shared pointer is not
  // owned anymore.
  if(self->structure) self->structure.reset(); 
  return 0;
}

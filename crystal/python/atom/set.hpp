extern "C"
{
  //! Structure holding shared pointer to an atom.
  struct Set
  { 
    PyObject_HEAD
    PyObject* ptr_base;
    std::set< std::string > *ptr_set;
  };
  //! Returns Py_True if key is contained in the set.
  static int set_contains(Set* _in, PyObject* _key);
  //! Returns Py_True if nothing in _b is in _a.
  static PyObject* set_isdisjoint(Set* _a, PyObject* _b); 
  //! Returns Py_True if everything in _a is in _b.
  static PyObject* set_issubset(Set* _a, PyObject* _b);
  //! Returns Py_True if everything in _a is in _b and something in _b is not in _a.
  static PyObject* set_istruesubset(Set* _a, PyObject* _b);
  //! Returns Py_True if everything in _b is in _a.
  static PyObject* set_issuperset(Set* _a, PyObject* _b);
  //! Returns Py_True if everything in _b is in _a and something in _a is not in _b.
  static PyObject* set_istruesuperset(Set* _a, PyObject* _b);

  //! Converts python object to an std::set<std::string>.
  static PyObject* to_cpp_set_(PyObject* _a, std::set<std::string> &_out);
  //! Converts an  std::set<std::string> to a python set.
  static PyObject* from_cpp_set_(std::set<std::string> const &_set);
  //! Converts a Set to a python set.
  static PyObject* set_as_set(Set* _a) { return from_cpp_set_(*_a->ptr_set); }
  //! Returns union of two sets.
  static PyObject* set_union(Set* _a, PyObject* _b);
  //! Reflection of union_().
  static PyObject* set_update(Set* _a, PyObject* _b);
  //! Reflection of union_().
  static PyObject* set_ior(Set *_a, PyObject *_tuple);
  //! Returns intersection of two sets.
  static PyObject* set_intersection(Set* _a, PyObject* _b);
  //! Reflection of difference().
  static PyObject* set_difference_update(Set* _a, PyObject* _b);
  //! Reflection of difference().
  static PyObject* set_isub(Set* _a, PyObject* _b);
  //! Creates a copy of _a from which anything contained in _b has been removed.
  static PyObject* set_difference(Set* _a, PyObject* _b);
  //! Equivalent to std::set_symmetric_difference.
  static PyObject* set_symmetric_difference(Set* _a, PyObject* _b);
  //! Reflection of symmetric_difference.
  static PyObject* set_symmetric_difference_update(Set *_a, PyObject *_b);
  //! Reflection of symmetric_difference.
  static PyObject* set_ixor(Set *_a, PyObject *_b);
  //! Reflection of intersection().
  static PyObject* set_intersection_update(Set *_a, PyObject *_b);
  //! Reflection of intersection().
  static PyObject* set_iand(Set *_a, PyObject *_b);
  //! Adds key to set.
  static PyObject* set_add(Set* _a, PyObject* _key);
  //! Removes key from set.
  static PyObject* set_remove(Set *_a, PyObject* _key);
  //! Removes and returns key from set.
  static PyObject* set_pop(Set *_a, PyObject* _key);
  //! Clears occupation set.
  static PyObject* set_clear(Set *_a) { _a->ptr_set->clear(); }
  //! \brief Discards key from set.
  //! \details Unlike ``pop`` or ``remove``, does not raise python exception if
  //!          key does not exist.
  static PyObject* set_discard(Set *_a, PyObject* _key);
  //! Python representation of set.
  static PyObject* set_repr(Set* _a);
  //! String representation of set.
  static PyObject* set_str(Set* _a);
  //! Returns Py_True if ``_b`` == ``_a``.
  static PyObject* set_equal(Set* _a, PyObject *_b);
  //! Returns Py_False if ``_b`` == ``_a``.
  static PyObject* set_nequal(Set* _a, PyObject *_b);
  //! Function to deallocate a string atom.
  static void set_dealloc(Set *_self);
  //! Returns size of set.
  static Py_ssize_t set_size(Set* _in);

  //! Rich comparison for sets.
  static PyObject *set_richcmp(PyObject *_a, PyObject *_b, int _op);

  //! \brief Creates new set. Debugging only.
  //! \details Sets should not be created from scratch. They are meant to wrap
  //!          atomic site occupations only.
  static PyObject* set_new_set(PyObject*, PyObject*);
  //! \brief Returns union of two sets.
  //! \details Checks which a and b is a Set from this extension.
  static PyObject* set_or(PyObject* _a, PyObject* _b);
  //! \brief Returns intersection of two sets.
  //! \details Checks which a and b is a Set from this extension.
  static PyObject* set_and(PyObject* _a, PyObject* _b);
  //! \brief Returns symmetric difference of two sets.
  //! \details Checks which a and b is a Set from this extension.
  static PyObject* set_xor(PyObject* _a, PyObject* _b);
  //! \brief Returns difference of two sets.
  //! \details Checks which a and b is a Set from this extension.
  static PyObject* set_sub(PyObject* _a, PyObject* _b);
};

#include "set_iterator.hpp"

static int set_contains(Set* _in, PyObject* _key)
{
  char *const string = PyString_AsString(_key);
  if(string == NULL) return -1;
  if( _in->ptr_set->count(string) == size_t(1) ) return 1;
  return 0;
}

static PyObject* set_isdisjoint(Set* _a, PyObject* _b)
{
  PyObject* iterator = PyObject_GetIter(_b);
  if(iterator == NULL) return NULL;
  while(PyObject* item = PyIter_Next(iterator))
  {
    char * const string = PyString_AsString(item);
    Py_DECREF(item);
    if(string == NULL) { Py_DECREF(iterator); return NULL; }
    if(_a->ptr_set->count(PyString_AS_STRING(item)) == size_t(1)) { Py_RETURN_FALSE; }
  }
  Py_DECREF(iterator);
  Py_RETURN_TRUE;
}

static PyObject* set_issubset(Set* _a, PyObject* _b)
{
  std::set<std::string>::const_iterator i_first = _a->ptr_set->begin();
  std::set<std::string>::const_iterator const i_end = _a->ptr_set->end();
  char method[] = "__contains__";
  char format[] = "s";
  for(; i_first != i_end; ++i_first)
  {
    PyObject *result = PyObject_CallMethod(_b, method, format, i_first->c_str());
    if(result == NULL) return NULL;
    Py_DECREF(result);
    if(result == Py_False) { Py_RETURN_FALSE; }
  }
  Py_RETURN_TRUE;
}
static PyObject* set_istruesubset(Set* _a, PyObject* _b)
{
  if( Py_ssize_t(_a->ptr_set->size()) >= PySequence_Size(_b) ) Py_RETURN_FALSE;
  return set_issubset(_a, _b);
}
static PyObject* set_issuperset(Set* _a, PyObject* _b)
{
  PyObject* iterator = PyObject_GetIter(_b);
  if(iterator == NULL) return NULL; 
  while(PyObject* item = PyIter_Next(iterator))
  {
    char * const string = PyString_AsString(item);
    Py_DECREF(item);
    if(string == NULL) {Py_DECREF(iterator); return NULL; }
    if(_a->ptr_set->count(PyString_AS_STRING(item)) == size_t(0)) { Py_RETURN_FALSE; }
  }
  Py_DECREF(iterator);
  Py_RETURN_TRUE;
}
static PyObject* set_istruesuperset(Set* _a, PyObject* _b)
{
  if( Py_ssize_t(_a->ptr_set->size()) <= PySequence_Size(_b) ) Py_RETURN_FALSE;
  return set_issuperset(_a, _b);
}
static PyObject* set_equal(Set* _a, PyObject *_b)
{
  if(set_issuperset(_a, _b) == Py_True and set_issubset(_a, _b) == Py_True){ Py_RETURN_TRUE; }
  Py_RETURN_FALSE;
}
static PyObject* set_nequal(Set* _a, PyObject *_b)
{
  if(set_issuperset(_a, _b) == Py_True and set_issubset(_a, _b) == Py_True){ Py_RETURN_FALSE; }
  Py_RETURN_TRUE;
}


static PyObject* to_cpp_set_(PyObject* _a, std::set<std::string> &_out)
{
  if(PyString_Check(_a))
  {
    _out.insert(PyString_AS_STRING(_a)); 
    Py_RETURN_NONE;
  }

  PyObject* iterator = PyObject_GetIter(_a);
  if(iterator == NULL) return NULL;
  while(PyObject* item = PyIter_Next(iterator))
  {
    char *const string = PyString_AsString(item);
    Py_DECREF(item);
    if(string == NULL) break;
    _out.insert(string);
  }
  Py_DECREF(iterator);
  if(PyErr_Occurred() != NULL) return NULL;
  Py_RETURN_NONE;
}
static PyObject* from_cpp_set_(std::set<std::string> const &_set)
{
  std::set<std::string>::const_iterator i_first = _set.begin();
  std::set<std::string>::const_iterator const i_end = _set.end();
  PyObject* result = PySet_New(NULL);
  if(result == NULL) return NULL;
  for(; i_first != i_end; ++i_first)
  {
    PyObject* string = PyString_FromString(i_first->c_str());
    if(string == NULL) { Py_DECREF(result); return NULL; }
    int const return_ = PySet_Add(result, string);
    Py_DECREF(string);
    if(return_ == 0) return NULL;
  }
  return result;
}

static PyObject* set_update(Set *_a, PyObject *_tuple)
{
  if(PyString_Check(_tuple) == 1) 
  {
    _a->ptr_set->insert(PyString_AS_STRING(_tuple));
    Py_RETURN_NONE;
  }
  PyObject* arg_iterator = PyObject_GetIter(_tuple);
  if(arg_iterator == NULL)
  { 
    std::cout << "error " << (PyErr_Occurred() != NULL) << "\n";
    return NULL;
  }
  while(PyObject* i_arg = PyIter_Next(arg_iterator))
  {
    char * const string = PyString_AsString(i_arg);
    Py_DECREF(i_arg);
    if(string == NULL) break; 
    _a->ptr_set->insert(string);
  }
  Py_DECREF(arg_iterator);
  if(PyErr_Occurred() != NULL) return NULL;
  Py_RETURN_NONE;
}
static PyObject* set_ior(Set *_a, PyObject *_tuple)
{
  PyObject * const result = set_update(_a, _tuple);
  if(result == NULL) return NULL;
  Py_XDECREF(result);
  Py_INCREF(_a);
  return (PyObject *)_a;
}

static PyObject* set_union(Set* _a, PyObject* _b)
{
  std::set<std::string> bset = *_a->ptr_set; 
  Set a; a.ptr_set = &bset;
  set_update(&a, _b);
  return from_cpp_set_(*a.ptr_set);
}

static PyObject* set_intersection_update(Set* _a, PyObject *_b)
{
  PyObject* arg_iterator = PyObject_GetIter(_b);
  if(arg_iterator == NULL) return NULL;
  while(PyObject* i_arg = PyIter_Next(arg_iterator))
  {
    if(PyString_Check(i_arg) == 1)
    {
      char * const key = PyString_AS_STRING(i_arg);
      bool const contains(_a->ptr_set->count(key) == 1);
      if(_a->ptr_set->size() > 1) _a->ptr_set->clear();
      if(contains) _a->ptr_set->insert(key);
    }
    else if(PyObject_HasAttrString(i_arg, "__contains__"))
    {
      std::set<std::string>::iterator i_first = _a->ptr_set->begin();
      std::set<std::string>::iterator const i_end = _a->ptr_set->end();
      char method[] = "__contains__";
      char format[] = "s";
      for(;i_first != i_end and PyErr_Occurred() == NULL; ++i_first)
      {
        PyObject * const result = PyObject_CallMethod(i_arg, method, format, i_first->c_str());
        if(result == NULL) break;
        if(result == Py_False) _a->ptr_set->erase(i_first);
        Py_DECREF(result);
      }
    }
    else LADA_PYERROR(TypeError, "Could not make sense of argument.");
    Py_DECREF(i_arg);
    if(PyErr_Occurred() != NULL) break;
  }
  Py_DECREF(arg_iterator);
  if(PyErr_Occurred() != NULL) return NULL;
  Py_RETURN_NONE;
}
static PyObject* set_iand(Set* _a, PyObject *_b)
{
  PyObject * const result = set_intersection_update(_a, _b);
  if(result == NULL) return NULL;
  Py_XDECREF(result);
  Py_INCREF(_a);
  return (PyObject *)_a;
}

static PyObject* set_intersection(Set* _a, PyObject* _b)
{
  std::set<std::string> bset = *_a->ptr_set; 
  Set a; a.ptr_set = &bset;
  PyObject *const result = set_intersection_update(&a, _b);
  if(result == NULL) return result;
  Py_DECREF(result);
  return from_cpp_set_(*a.ptr_set);
}

static PyObject* set_difference_update(Set *_a, PyObject *_b)
{
  PyObject* arg_iterator = PyObject_GetIter(_b);
  if(arg_iterator == NULL) return NULL;
  while(PyObject* i_arg = PyIter_Next(arg_iterator))
  {
    if(PyString_Check(i_arg) == 1)
    {
      std::set<std::string>::iterator const i_found 
        = _a->ptr_set->find(PyString_AS_STRING(i_arg));
      if(i_found != _a->ptr_set->end()) _a->ptr_set->erase(i_found);
    }
    else if(PyObject* seq_iterator = PyObject_GetIter(i_arg))
    {
      while(PyObject *i_str = PyIter_Next(seq_iterator))
      {
        char * const string = PyString_AsString(i_arg);
        Py_DECREF(i_str);
        if(string == NULL) break;
        std::set<std::string>::iterator const i_found 
          = _a->ptr_set->find(string);
        if(i_found != _a->ptr_set->end()) _a->ptr_set->erase(i_found);
      }
      Py_DECREF(seq_iterator);
    }
    else LADA_PYERROR(TypeError, "Could not make sense of argument."); 
    Py_DECREF(i_arg);
    if(PyErr_Occurred() != NULL) break;
  }
  Py_DECREF(arg_iterator);
  if(PyErr_Occurred() != NULL) return NULL;
  Py_RETURN_NONE;
}
static PyObject* set_isub(Set *_a, PyObject *_b)
{
  PyObject * const result = set_difference_update(_a, _b);
  if(result == NULL) return NULL;
  Py_DECREF(result);
  Py_INCREF(_a);
  return (PyObject *)_a;
}
static PyObject* set_difference(Set* _a, PyObject* _b)
{
  std::set<std::string> bset = *_a->ptr_set; 
  Set a; a.ptr_set = &bset;
  PyObject *const result = set_difference_update(&a, _b);
  if(result == NULL) return NULL;
  Py_DECREF(result);
  return from_cpp_set_(*a.ptr_set);
}

static PyObject* set_symmetric_difference_update(Set *_a, PyObject *_b)
{
  std::set<std::string> bset, result;
  if(PyObject* return_ = to_cpp_set_(_b, bset)) Py_DECREF(return_); 
  else return NULL;
  std::insert_iterator< std::set<std::string> > inserter(result, result.begin());
  std::set_symmetric_difference(_a->ptr_set->begin(), _a->ptr_set->end(), bset.begin(), bset.end(), inserter);
  *_a->ptr_set = result;
  Py_RETURN_NONE;
}
static PyObject* set_ixor(Set *_a, PyObject *_b)
{
  PyObject * const result = set_symmetric_difference_update(_a, _b);
  if(result == NULL) return NULL;
  Py_XDECREF(result);
  Py_INCREF(_a);
  return (PyObject *)_a;
}
static PyObject* set_symmetric_difference(Set* _a, PyObject* _b)
{
  std::set<std::string> aset;
  aset = *_a->ptr_set;
  Set result; result.ptr_set = &aset;
  PyObject* const return_ = set_ixor(&result, _b);
  if(return_ == NULL) return NULL;
  Py_DECREF(return_);
  return from_cpp_set_(aset);
}
static PyObject* set_add(Set* _a, PyObject* _key)
{
  if(char * const string = PyString_AsString(_key))
  {
    _a->ptr_set->insert(string); 
    Py_RETURN_NONE;
  }
  return NULL;
}
static PyObject* set_remove(Set *_a, PyObject* _key)
{
  if(char * const string = PyString_AsString(_key))
  {
    std::set<std::string>::iterator const i_found = _a->ptr_set->find(string);
    if(i_found == _a->ptr_set->end())
    {
      LADA_PYERROR(KeyError, "Atomic specie not in set.");
      return NULL;
    }
    else _a->ptr_set->erase(i_found);
    Py_RETURN_NONE;
  }
  return NULL;
}
static PyObject* set_pop(Set *_a, PyObject* _key)
{
  if(char * const string = PyString_AsString(_key))
  {
    std::set<std::string>::iterator const i_found = _a->ptr_set->find(string);
    if(i_found == _a->ptr_set->end())
      LADA_PYERROR(KeyError, "Atomic specie not in set.");
    else if(PyObject *result = PyString_FromString(i_found->c_str()))
    {
      _a->ptr_set->erase(i_found);
      return result;
    }
  }
  return NULL;
}
static PyObject* set_discard(Set *_a, PyObject* _key)
{
  if(char * const string = PyString_AsString(_key))
  {
    std::set<std::string>::iterator const i_found = _a->ptr_set->find(PyString_AS_STRING(_key));
    if(i_found != _a->ptr_set->end()) _a->ptr_set->erase(i_found);
    Py_RETURN_NONE;
  }
  return NULL;
}

static PyObject* set_repr(Set* _a)
{
  if( _a->ptr_set->size() == 0 ) return PyString_FromString("set()");
  std::ostringstream sstr;
  std::set<std::string>::const_iterator i_first = _a->ptr_set->begin();
  std::set<std::string>::const_iterator const i_end = _a->ptr_set->end();
  sstr << "set([\"" << *i_first << "\"";
  for(++i_first; i_first != i_end; ++i_first)
    sstr << ", \"" << *i_first << "\"";
  sstr << "])";
  return PyString_FromString(sstr.str().c_str());
}
static PyObject* set_str(Set* _a)
{
  if( _a->ptr_set->size() == 0 ) return PyString_FromString("set()");
  std::ostringstream sstr;
  std::set<std::string>::const_iterator i_first = _a->ptr_set->begin();
  std::set<std::string>::const_iterator const i_end = _a->ptr_set->end();
  sstr << "{\"" << *i_first << "\"";
  for(++i_first; i_first != i_end; ++i_first)
    sstr << ", \"" << *i_first << "\"";
  sstr << "}";
  return PyString_FromString(sstr.str().c_str());
}
// Function to deallocate a string atom.
static void set_dealloc(Set *_self)
{
  if(_self->ptr_base != NULL)
  {
    PyObject *dummy = _self->ptr_base;
    _self->ptr_base = NULL;
    Py_DECREF(dummy);
    _self->ptr_set = NULL;
  }
  else if(_self->ptr_set != NULL) { delete _self->ptr_set; }
  _self->ob_type->tp_free((PyObject*)_self);
}

static Py_ssize_t set_size(Set* _in)
  { return Py_ssize_t(_in->ptr_set->size()); }

static PyMethodDef set_methods[] = {
    {"isdisjoint", (PyCFunction)set_isdisjoint, METH_O,
     "Returns true if self and input set have nothing in common." },
    {"issubset", (PyCFunction)set_issubset, METH_O,
     "Returns true if self is a subset of the input set." },
    {"issuperset", (PyCFunction)set_issuperset, METH_O,
     "Returns true if self is a superset of the input set." },
    {"union", (PyCFunction)set_union, METH_O,
     "Returns union of self and input set." },
    {"intersection", (PyCFunction)set_intersection, METH_VARARGS,
     "Returns intersection of self and input set." },
    {"difference", (PyCFunction)set_difference, METH_VARARGS,
     "Returns difference of self and input set." },
    {"symmetric_difference", (PyCFunction)set_symmetric_difference, METH_O,
     "Returns symmetric_difference of self and input set." },
    {"copy", (PyCFunction)set_as_set, METH_NOARGS, 
     "Returns a copy of self (as python intrinsic type)." },
    {"update", (PyCFunction)set_update, METH_O,
     "Adds objects in sequence to self." },
    {"difference_update", (PyCFunction)set_difference_update, METH_VARARGS,
     "Removes objects in input from self." },
    {"intersection_update", (PyCFunction)set_intersection_update, METH_VARARGS,
     "self becomes intersection of self and input sequences." },
    {"symmetric_difference_update", (PyCFunction)set_symmetric_difference_update, METH_VARARGS,
     "self becomes symmetric difference of self and input sequences." },
    {"add", (PyCFunction)set_add, METH_O, "Adds key to occupation set." },
    {"remove", (PyCFunction)set_remove, METH_O, "Removes key to occupation set." },
    {"discard", (PyCFunction)set_discard, METH_O,
      "Removes key to occupation set.\n\n"
      "Unlike ``pop`` or ``remove``, does not throw if occupation is not present." },
    {"pop", (PyCFunction)set_pop, METH_O, "Removes key to occupation set and returns it." },
    {"clear", (PyCFunction)set_clear, METH_NOARGS, "Clears occupation set." },
    {NULL}  /* Sentinel */
};

static PySequenceMethods set_as_sequence = {
    (lenfunc) set_size, 0, 0, 0, 0, 0, 0,
    (objobjproc) set_contains, 0, 0
  };

static PyNumberMethods set_as_number = { 
     (binaryfunc) set_or, (binaryfunc) set_sub,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
     binaryfunc (set_and),
     binaryfunc (set_xor),
     binaryfunc (set_or),
     0, 0, 0, 0, 0, 0, 
     binaryfunc (set_ior),
     binaryfunc (set_isub),
     0, 0, 0, 0, 0, 0, 
     binaryfunc (set_iand),
     binaryfunc (set_ixor),
     binaryfunc (set_ior),
     0, 0, 0, 0, 0
}; 

  

static PyTypeObject set_type = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "atom.Set",                /*tp_name*/
    sizeof(Set),               /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)set_dealloc,   /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    (reprfunc)set_repr,        /*tp_repr*/
    &set_as_number,
    &set_as_sequence,          /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    (reprfunc)set_str,         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_CHECKTYPES,
    "Atomic occupation of a site.\n\n"
    "Wrapper around a C++ object. It works exactly like a python set, with the "
    "restriction that it can only contain strings. It cannot be created. It always "
    "references the occupation of an atom.",
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    (richcmpfunc)set_richcmp,  /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    (getiterfunc)setiterator_create,  /* tp_iter */
    0,		               /* tp_iternext */
    set_methods                /* tp_methods */
};

static PyObject* set_new_set(PyObject* _module, PyObject* _in)
{
  Set* result = NULL;
  PyObject *input = NULL;
  if(not PyArg_UnpackTuple(_in, "_new_set", 0, 1, &input)) goto error;
  // creates new object and check result.
  result = PyObject_New(Set, (&set_type));
  if(result == NULL) goto error;
  
  // on success initializes object.
  result->ptr_base = NULL;
  result->ptr_set = new (std::nothrow) std::set<std::string>;
  if(result->ptr_set == NULL) 
  {
    LADA_PYERROR(internal, "Could not create set.");
    goto error;
  }
  if(input != NULL)
  {
    if(PyObject* return_ = to_cpp_set_(input, *result->ptr_set)) Py_DECREF(return_);
    else goto error;
  }
  
  return (PyObject*)result;

  error:
    Py_XDECREF(result);
    return NULL;
};

static PyObject* set_or(PyObject* _a, PyObject* _b)
{
  if(_a->ob_type == &set_type) { return set_union((Set*)_a, _b); }
  if(_b->ob_type == &set_type) { return set_union((Set*)_b, _a); }
  LADA_PYERROR(TypeError, "Neither input is a lada.crystal.atom.Set");
  return NULL;
}
static PyObject* set_and(PyObject* _a, PyObject* _b)
{
  if(_a->ob_type == &set_type) { return set_intersection((Set*)_a, _b); }
  if(_b->ob_type == &set_type) { return set_intersection((Set*)_b, _a); }
  LADA_PYERROR(TypeError, "Neither input is a lada.crystal.atom.Set");
  return NULL;
}
static PyObject* set_xor(PyObject* _a, PyObject* _b)
{
  if(_a->ob_type == &set_type) { return set_symmetric_difference((Set*)_a, _b); }
  if(_b->ob_type == &set_type) { return set_symmetric_difference((Set*)_b, _a); }
  LADA_PYERROR(TypeError, "Neither input is a lada.crystal.atom.Set");
  return NULL;
}
static PyObject* set_sub(PyObject* _a, PyObject* _b)
{
  if(_a->ob_type == &set_type) { return set_difference((Set*)_a, _b); }
  if(_b->ob_type != &set_type)
  {
    LADA_PYERROR(TypeError, "Neither input is a lada.crystal.atom.Set");
    return NULL;
  }
  // set pointer to C++ set.
  std::set<std::string> const * const b = ((Set*)_b)->ptr_set;
  // construct result set.
  PyObject* result = PySet_New(NULL);
  if(result == NULL) return NULL;
  // loop over objects in a. Adds to result if not in b.
  PyObject* arg_iterator = PyObject_GetIter(_a);
  if(arg_iterator == NULL) { Py_DECREF(result); return NULL; }
  char method[] = "add";
  char format[] = "O";
  while(PyObject* i_arg = PyIter_Next(arg_iterator))
  {
    char * const string = PyString_AsString(i_arg);
    if(string == NULL) { Py_DECREF(i_arg); break; }
    if(b->count(string) == 0) Py_XDECREF(PyObject_CallMethod(result, method, format, i_arg));
    Py_DECREF(i_arg);
    if(PyErr_Occurred() != NULL) break;
  }
  Py_DECREF(arg_iterator);
  if(PyErr_Occurred() == NULL) return result;

  Py_DECREF(result);
  return NULL;
}

static PyObject *set_richcmp(PyObject *_a, PyObject *_b, int _op)
{
  PyObject *result = NULL;
  switch(_op)
  {
    case Py_LT: result = set_istruesubset((Set*)_a, _b); break;
    case Py_LE: result = set_issubset((Set*)_a, _b); break;
    case Py_EQ: result = set_equal((Set*)_a, _b); break;
    case Py_NE: result = set_nequal((Set*)_a, _b); break;
    case Py_GT: result = set_istruesuperset((Set*)_a, _b); break;
    case Py_GE: result = set_issuperset((Set*)_a, _b); break;
    default: break;
  };
  if(result == NULL) LADA_PYERROR(internal, "Unknown rich comparison.");
  return result;
}


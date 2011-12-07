extern "C"
{
  //! Returns Py_True if key is contained in the set.
  static int sequence_contains(Sequence* _in, PyObject* _key);
  //! Converts python object to an std::vector<std::string>.
  static PyObject* to_cpp_sequence_(PyObject* _a, std::vector<std::string> &_out);
  //! Converts an  std::vector<std::string> to a python set.
  static PyObject* from_cpp_sequence_(std::vector<std::string> const &_seq);
  //! Converts a Sequence to a python list.
  static PyObject* sequence_as_list(Sequence* _self) { return from_cpp_sequence_(*_self->ptr_seq); }
  //! Adds key to set.
  static PyObject* sequence_append(Sequence* _self, PyObject* _key);
  //! Removes key from set.
  static PyObject* sequence_remove(Sequence *_self, PyObject* _key);
  //! Removes and returns key from set.
  static PyObject* sequence_pop(Sequence *_self, PyObject* _key);
  //! Clears occupation set.
  static PyObject* sequence_clear(Sequence *_self) { _self->ptr_seq->clear(); Py_RETURN_NONE;}
  //! Python representation of set.
  static PyObject* sequence_repr(Sequence* _self);
  //! Function to deallocate a string atom.
  static void sequence_dealloc(Sequence *_self);
  //! Returns size of set.
  static Py_ssize_t sequence_size(Sequence* _in) { return Py_ssize_t(_in->ptr_seq->size()); }
  //! Concatenate two sequences.
  static PyObject* sequence_concat(Sequence* _self, PyObject *_b);
  //! Repeats a sequence.
  static PyObject* sequence_repeat(Sequence* _self, Py_ssize_t _n);
  //! In-place concatenation.
  static PyObject* sequence_inplaceconcat(Sequence* _self, PyObject *_b);
  //! In-place repeat.
  static PyObject* sequence_inplacerepeat(Sequence* _self, Py_ssize_t _n);
  //! Retrives item. No slice.
  static PyObject* sequence_getitem(Sequence *_self, Py_ssize_t _index);
  //! Sets/deletes item. No slice.
  static int sequence_setitem(Sequence *_self, Py_ssize_t _index, PyObject *_replace);

  //! \brief Creates new set. Debugging only.
  //! \details Sets should not be created from scratch. They are meant to wrap
  //!          atomic site occupations only.
  static PyObject* sequence_new_set(PyObject*, PyObject*);
  // Helper function for equality.
  static int sequence_equal_impl(Sequence* _self, PyObject *_b);
  //! Rich (lexicographical) comparison for sequences.
  static PyObject *sequence_richcmp(PyObject *_self, PyObject *_b, int _op);
  //! Appends specie to list.
  static PyObject* sequence_append(Sequence* _self, PyObject *_b);
  //! Extends species with other list of species.
  static PyObject* sequence_extend(Sequence* _self, PyObject *_b);
  //! Insert specie at given position
  static PyObject* sequence_insert(Sequence* _self, PyObject *_b);
  //! Returns index of first matching specie.
  static PyObject* sequence_index(Sequence* _self, PyObject *_b);
  //! Returns number of matching species in list.
  static PyObject* sequence_count(Sequence* _self, PyObject *_b);
  //! Sorts list in place.
  static PyObject* sequence_sort(Sequence* _self);
  //! Reverse list in place.
  static PyObject* sequence_reverse(Sequence* _self);
  //! garbage collection.
  static int sequence_traverse(Sequence *self, visitproc visit, void *arg);
  //! garbage collection.
  static int sequence_gcclear(Sequence *_self);
};

#include "sequence_iterator.hpp"

static int sequence_contains(Sequence* _in, PyObject* _key)
{
  char *const string = PyString_AsString(_key);
  if(string == NULL) return -1;
  if(std::find(_in->ptr_seq->begin(), _in->ptr_seq->end(), string) != _in->ptr_seq->end()) return 1;
  return 0;
}

static PyObject* from_cpp_sequence_(std::vector<std::string> const &_seq)
{
  std::vector<std::string>::const_iterator i_first = _seq.begin();
  std::vector<std::string>::const_iterator const i_end = _seq.end();
  PyObject* result = PyList_New(_seq.size());
  if(result == NULL) return NULL;
  for(Py_ssize_t i(0); i_first != i_end; ++i_first, ++i)
  {
    PyObject* string = PyString_FromString(i_first->c_str());
    if(string == NULL) { Py_DECREF(result); return NULL; }
    PyList_SET_ITEM(result, i, string);
  }
  return result;
}

static PyObject* set_append(Sequence* _self, PyObject* _key)
{
  if(char * const string = PyString_AsString(_key))
  {
    _self->ptr_seq->push_back(string); 
    Py_RETURN_NONE;
  }
  return NULL;
}

static inline long get_index(Sequence *_self, PyObject* _index)
{
  long index = PyInt_AsLong(_index);
  if(index == -1 and PyErr_Occurred() != NULL) return -1;
  if(index < 0) index += _self->ptr_seq->size();
  if(index < 0 or index >= _self->ptr_seq->size())
  {
    LADA_PYERROR(IndexError, "Index out of range.");
    return -1;
  }
  return index;
}
static PyObject* sequence_pop(Sequence *_self, PyObject* _index)
{
  long const index = get_index(_self, _index);
  if(index < 0) return NULL;
  if(index + 1 == _self->ptr_seq->size())
  {
    std::string const string = _self->ptr_seq->back();
    _self->ptr_seq->pop_back();
    return PyString_FromString(string.c_str());
  }
  PyObject *result = PyString_FromString(_self->ptr_seq->operator[](index).c_str());
  _self->ptr_seq->erase(_self->ptr_seq->begin() + index);
  return result;
}
static int sequence_setitem(Sequence *_self, Py_ssize_t _index, PyObject* _replace)
{
  long const index = _index < 0 ? _index + _self->ptr_seq->size(): _index;
  if(index < 0 or index >= _self->ptr_seq->size())
  { 
    LADA_PYERROR(IndexError, "Index is out of range in item assignement.");
    return -1;
  }
  if(_replace == NULL)
  {
    if(index+1 == _self->ptr_seq->size()) _self->ptr_seq->pop_back();
    else _self->ptr_seq->erase(_self->ptr_seq->erase(_self->ptr_seq->begin()+index));
  }
  else 
  {
    char *const string = PyString_AsString(_replace);
    if(string == NULL) return -1;
    (*_self->ptr_seq)[index] = string;
  }
  return 0;
}
static PyObject* sequence_getitem(Sequence *_self, Py_ssize_t _index)
{
  long const index = _index < 0 ? _index + _self->ptr_seq->size(): _index;
  if(index < 0 or index >= _self->ptr_seq->size())
  { 
    LADA_PYERROR(IndexError, "Index is out of range.");
    return NULL;
  }
  return PyString_FromString(_self->ptr_seq->operator[](index).c_str());
}

static PyObject* sequence_repr(Sequence* _self)
{
  if( _self->ptr_seq->size() == 0 ) return PyString_FromString("[]");
  std::ostringstream sstr;
  std::vector<std::string>::const_iterator i_first = _self->ptr_seq->begin();
  std::vector<std::string>::const_iterator const i_end = _self->ptr_seq->end();
  sstr << "['" << *i_first << "'";
  for(++i_first; i_first != i_end; ++i_first)
    sstr << ", '" << *i_first << "'";
  sstr << "]";
  return PyString_FromString(sstr.str().c_str());
}
// Function to deallocate a string atom.
static void sequence_dealloc(Sequence *_self)
{
  _self->ptr_atom.reset();
  _self->ptr_seq = NULL;
  _self->ob_type->tp_free((PyObject*)_self);
}

// Concatenate two sequences.
static PyObject* sequence_concat(Sequence* _self, PyObject *_b)
{
  PyObject *result = from_cpp_sequence_(*_self->ptr_seq);
  if(result == NULL) return NULL;

  PyObject* iterator = PyObject_GetIter(_b);
  if(iterator == NULL) {Py_DECREF(result); return NULL;}
  while(PyObject* item = PyIter_Next(iterator))
  {
    bool const error = PyList_Append(result, item) == -1;
    Py_DECREF(item);
    if(error) {Py_DECREF(iterator); Py_DECREF(result); return NULL;}
  }
  Py_DECREF(iterator);
  return result;
}
static PyObject* sequence_repeat(Sequence* _self, Py_ssize_t _n)
{
  Py_ssize_t const N(_self->ptr_seq->size());
  PyObject* result = PyList_New(N * _n);
  if(result == NULL) return NULL;
  if(_n == 0 or N == 0) return result;
  std::vector<std::string>::const_iterator const i_end = _self->ptr_seq->end();
  std::vector<std::string>::const_iterator i_first = _self->ptr_seq->begin();
  for(Py_ssize_t j(0); i_first != i_end; ++i_first, ++j)
  {
    PyObject* string = PyString_FromString(i_first->c_str());
    if(string == NULL) { Py_DECREF(result); return NULL; }
    PyList_SET_ITEM(result, j, string);
    for(Py_ssize_t i(1); i < _n; ++i) 
    {
      Py_INCREF(string);
      PyList_SET_ITEM(result, j + i * N, string);
    }
  }
  return result;
}
static PyObject* sequence_inplaceconcat(Sequence* _self, PyObject* _b)
{
  std::back_insert_iterator< std::vector<std::string> > i_inserter(*_self->ptr_seq);
  PyObject* iterator = PyObject_GetIter(_b);
  if(iterator == NULL) {Py_DECREF(iterator); return NULL;}
  while(PyObject* item = PyIter_Next(iterator))
  {
    char *const string = PyString_AsString(item);
    Py_DECREF(item);
    if(string == 0) { Py_DECREF(iterator); return NULL; }
    *i_inserter = string;
    ++i_inserter;
  }
  Py_DECREF(iterator);
  Py_INCREF(_self);
  return (PyObject*)_self;
}
static PyObject* sequence_inplacerepeat(Sequence* _self, Py_ssize_t _n)
{
  std::vector< std::string> :: size_type N(_self->ptr_seq->size());
  if(_n > 1 and N != 0)
  {
    _self->ptr_seq->reserve(N * _n);
    std::back_insert_iterator< std::vector<std::string> > i_inserter(*_self->ptr_seq);
    std::vector<std::string>::const_iterator const i_first = _self->ptr_seq->begin();
    std::vector<std::string>::const_iterator const i_end = _self->ptr_seq->end();
    for(Py_ssize_t i(1); i < _n; ++i) std::copy(i_first, i_end, i_inserter);
  }
  Py_INCREF(_self);
  return (PyObject*)_self;
}

static PyObject* sequence_append(Sequence* _self, PyObject *_b)
{
  char * const string = PyString_AsString(_b);
  if(string == NULL) return NULL;
  _self->ptr_seq->push_back(string);
  Py_RETURN_NONE;
}
static PyObject* sequence_insert(Sequence* _self, PyObject *_b)
{
  char *specie;
  int index;
  if(not PyArg_ParseTuple(_b, "is", &index, &specie) ) return NULL;
  if(index < 0) index += _self->ptr_seq->size();
  if(index < 0 or index > _self->ptr_seq->size())
  {
    LADA_PYERROR(IndexError, "Insertion index is out of range.");
    return NULL;
  }
  _self->ptr_seq->insert(_self->ptr_seq->begin() + index, specie);
  Py_RETURN_NONE;
}
static PyObject* sequence_index(Sequence* _self, PyObject *_b)
{
  char * const string = PyString_AsString(_b);
  if(string == NULL) return NULL;
  std::vector<std::string> :: const_iterator const i_found 
     = std::find(_self->ptr_seq->begin(), _self->ptr_seq->end(), string);
  if(i_found == _self->ptr_seq->end()) 
  {
    LADA_PYERROR(ValueError, (std::string(string) + " is not in specie list.").c_str());
    return NULL;
  }
  long const result = i_found - _self->ptr_seq->begin(); 
  return PyInt_FromLong(result);
}
static PyObject* sequence_count(Sequence* _self, PyObject *_b)
{
  char * const string = PyString_AsString(_b);
  if(string == NULL) return NULL;
  long const result = std::count(_self->ptr_seq->begin(), _self->ptr_seq->begin(), string);
  return PyInt_FromLong(result);
}
static PyObject* sequence_sort(Sequence* _self)
{
  std::sort(_self->ptr_seq->begin(), _self->ptr_seq->end());
  Py_RETURN_NONE;
}
static PyObject* sequence_reverse(Sequence* _self)
{
  std::reverse(_self->ptr_seq->begin(), _self->ptr_seq->end());
  Py_RETURN_NONE;
}


static PyMethodDef sequence_methods[] = {
    {"copy", (PyCFunction)sequence_as_list, METH_NOARGS, 
     "Returns a copy of the sequence as a python list." },
    {"append", (PyCFunction)sequence_append, METH_O, "Appends specie to list." },
    {"extend", (PyCFunction)sequence_extend, METH_O, "Extend list with another list of species." },
    {"insert", (PyCFunction)sequence_insert, METH_VARARGS, "Insert specie at given position." },
    {"pop", (PyCFunction)sequence_pop, METH_O|METH_COEXIST, "Removes key to occupation set and returns it." },
    {"index", (PyCFunction)sequence_index, METH_O, "Returns index of first matching specie." },
    {"count", (PyCFunction)sequence_count, METH_NOARGS,  "Counts instances of a specie in a list." }, 
    {"reverse", (PyCFunction)sequence_reverse, METH_NOARGS, "Reverse list of species in place." }, 
    {"sort", (PyCFunction)sequence_sort, METH_NOARGS, "Sorts (lexicographically) list of species in place." }, 
    {"clear", (PyCFunction)sequence_clear, METH_NOARGS, "Clears occupation set." },
    {NULL}  /* Sentinel */
};

static PySequenceMethods sequence_as_sequence = {
    (lenfunc) sequence_size,
    (binaryfunc) sequence_concat,
    (ssizeargfunc) sequence_repeat,
    (ssizeargfunc) sequence_getitem,
    0,
    (ssizeobjargproc) sequence_setitem,
    0,
    (objobjproc) sequence_contains,
    (binaryfunc) sequence_inplaceconcat,
    (ssizeargfunc) sequence_inplacerepeat
};

PyTypeObject sequence_type = {
    PyObject_HEAD_INIT(NULL)
    0,                                                       /*ob_size*/
    "atom.Sequence",                                         /*tp_name*/
    sizeof(Sequence),                                        /*tp_basicsize*/
    0,                                                       /*tp_itemsize*/
    (destructor)sequence_dealloc,                            /*tp_dealloc*/
    0,                                                       /*tp_print*/
    0,                                                       /*tp_getattr*/
    0,                                                       /*tp_setattr*/
    0,                                                       /*tp_compare*/
    (reprfunc)sequence_repr,                                 /*tp_repr*/
    0,
    &sequence_as_sequence,                                   /*tp_as_sequence*/
    0,                                                       /*tp_as_mapping*/
    0,                                                       /*tp_hash */
    0,                                                       /*tp_call*/
    (reprfunc)sequence_repr,                                 /*tp_str*/
    0,                                                       /*tp_getattro*/
    0,                                                       /*tp_setattro*/
    0,                                                       /*tp_as_buffer*/
    Py_TPFLAGS_HAVE_SEQUENCE_IN | Py_TPFLAGS_HAVE_INPLACEOPS | Py_TPFLAGS_HAVE_RICHCOMPARE
      | Py_TPFLAGS_HAVE_ITER | Py_TPFLAGS_HAVE_CLASS,
    "Atomic occupation of a site.\n\n"
    "Wrapper around a C++ object. It works exactly like a python list, with the "
    "restriction that it can only contain strings. It cannot be created. It always "
    "references the occupation of an atom.",
    0,//(traverseproc)sequence_traverse,                         /* tp_traverse */
    0,//(inquiry)sequence_gcclear,                               /* tp_clear */
    (richcmpfunc)sequence_richcmp,                           /* tp_richcompare */
    0,		                                             /* tp_weaklistoffset */
    (getiterfunc)sequenceiterator_create,                    /* tp_iter */
    0,		                                             /* tp_iternext */
    sequence_methods                                         /* tp_methods */
};

static PyObject* sequence_new_sequence(PyObject* _module, PyObject* _in)
{
  Sequence* result = NULL;
  PyObject *input = NULL;
  if(not PyArg_UnpackTuple(_in, "_new_sequence", 0, 1, &input)) goto error;
  // creates new object and check result.
  result = PyObject_New(Sequence, &sequence_type);
  if(result == NULL) goto error;
  
  typedef LaDa::crystal::AtomData< std::vector<std::string> > t_Atom;
  new(&result->ptr_atom) boost::shared_ptr<t_Atom>;
  // on success initializes object.
  result->ptr_seq = new (std::nothrow) std::vector<std::string>;
  if(result->ptr_seq == NULL) 
  {
    LADA_PYERROR(internal, "Could not create sequence.");
    goto error;
  }
  if(input != NULL)
  {
    if(PyObject* return_ = to_cpp_sequence_(input, *result->ptr_seq)) Py_DECREF(return_);
    else goto error;
  }
  
  return (PyObject*)result;

  error:
    Py_XDECREF(result);
    return NULL;
};

static inline int sequence_cmp_object(Sequence *_self, PyObject *_b)
{
  if(PyString_Check(_b) == 1)
  {
    char const * const string = PyString_AsString(_b);
    if(_self->ptr_seq->size() == 0) return std::string(string) == "" ? 0: -1;
    if(_self->ptr_seq->front() == string) return 0;
    return _self->ptr_seq->front() < string ? -1: 1;
  }
  PyObject* iterator = PyObject_GetIter(_b);
  if(iterator == NULL) return -2;
  // loop until we find a strictly < or strictly > comparison.
  std::vector<std::string>::const_iterator i_first = _self->ptr_seq->begin();
  std::vector<std::string>::const_iterator const i_end = _self->ptr_seq->end();
  bool empty = true; 
  while(PyObject* item = PyIter_Next(iterator))
  {
    empty = false;
    if(i_first == i_end) { Py_DECREF(item); Py_DECREF(iterator); return -1; }
    char *const string = PyString_AsString(item);
    Py_DECREF(item);
    if(string == NULL) { Py_DECREF(iterator); return -2; }
    if(*i_first == string) { ++i_first; continue; }

    // found true. comparison 
    Py_DECREF(iterator);
    return *i_first < string ? -1: 1;
  }
  Py_DECREF(iterator);
  if(empty) return i_first == i_end ? 0: 1;
  return i_first != i_end ? 1: 0;
}
static inline int sequence_cmp_sequence( std::vector<std::string> const &_self,
                                         std::vector<std::string> const &_b )
{
  size_t const Na(_self.size());
  size_t const Nb(_b.size());
  if(Na == Nb) 
  {
    typedef std::vector<std::string> :: const_iterator t_cit;
    typedef std::pair<t_cit, t_cit> t_pair;
    std::vector<std::string> :: const_iterator const i_a_end = _self.end();
    std::vector<std::string> :: const_iterator const i_b_end = _b.end();
    t_pair const result = std::mismatch(_self.begin(), _self.end(), _b.begin());
    if(result.first == i_a_end) return 0;
    return std::lexicographical_compare(result.first, i_a_end, result.second, i_b_end) ? -1: 1;
  }
  return std::lexicographical_compare(_self.begin(), _self.end(), _b.begin(), _b.end()) ? -1: 1;
}

static PyObject *sequence_richcmp(PyObject *_self, PyObject *_b, int _op)
{
  // swap operations and object if a is not a sequence.
  if(_self->ob_type != &sequence_type)
  { 
    if(_b->ob_type != &sequence_type) 
      LADA_PYERROR(internal, "Neither object is a sequence. Shouldn't be here.");
    switch(_op)
    {
      case Py_LT: _op = Py_GE; break;
      case Py_LE: _op = Py_GT; break;
      case Py_EQ: _op = Py_EQ; break;
      case Py_NE: _op = Py_NE; break;
      case Py_GT: _op = Py_LE; break;
      case Py_GE: _op = Py_LT; break;
      default: break;
    };
    std::swap(_self, _b); 
  }
  size_t const return_ = _b->ob_type == &sequence_type ?
    sequence_cmp_sequence(*((Sequence*)_self)->ptr_seq, *((Sequence*)_b)->ptr_seq):
    sequence_cmp_object((Sequence*)_self, _b);
  if(return_ == -2) return NULL;
  PyObject* result = NULL;
  switch(_op)
  {
    case Py_LT: result = return_ == -1 ? Py_True: Py_False; Py_INCREF(result); break;
    case Py_LE: result = return_ !=  1 ? Py_True: Py_False; Py_INCREF(result); break;
    case Py_EQ: result = return_ ==  0 ? Py_True: Py_False; Py_INCREF(result); break;
    case Py_NE: result = return_ !=  0 ? Py_True: Py_False; Py_INCREF(result); break;
    case Py_GT: result = return_ ==  1 ? Py_True: Py_False; Py_INCREF(result); break;
    case Py_GE: result = return_ != -1 ? Py_True: Py_False; Py_INCREF(result); break;
    default: LADA_PYERROR(internal, "Unknown rich comparison."); break;
  };
  return result;
}

static PyObject* sequence_extend(Sequence* _self, PyObject* _b)
{
  std::back_insert_iterator< std::vector<std::string> > i_inserter(*_self->ptr_seq);
  if(_b->ob_type == &sequence_type)
  {
    std::vector<std::string> const &b = *((Sequence*)_b)->ptr_seq;
    _self->ptr_seq->reserve(_self->ptr_seq->size() + b.size());
    std::copy(b.begin(), b.end(), std::back_inserter(*_self->ptr_seq));
  }
  else
  {
    PyObject* iterator = PyObject_GetIter(_b);
    if(iterator == NULL) {Py_DECREF(iterator); return NULL;}
    while(PyObject* item = PyIter_Next(iterator))
    {
      char *const string = PyString_AsString(item);
      Py_DECREF(item);
      if(string == 0) { Py_DECREF(iterator); return NULL; }
      *i_inserter = string;
      ++i_inserter;
    }
    Py_DECREF(iterator);
  }
  Py_RETURN_NONE;
}

static PyObject* to_cpp_sequence_(PyObject* _a, std::vector<std::string> &_out)
{
  if(PyString_Check(_a))
  {
    _out.push_back(PyString_AS_STRING(_a));
    Py_RETURN_NONE;
  }
  else if(_a->ob_type == &sequence_type)
  {
    _out = *((Sequence*)_a)->ptr_seq;
    Py_RETURN_NONE;
  }
  PyObject* iterator = PyObject_GetIter(_a);
  if(iterator == NULL) return NULL;
  while(PyObject* item = PyIter_Next(iterator))
  {
    char *const string = PyString_AsString(item);
    Py_DECREF(item);
    if(string == NULL) break;
    _out.push_back(string);
  }
  Py_DECREF(iterator);
  if(PyErr_Occurred() != NULL) return NULL;
  Py_RETURN_NONE;
}

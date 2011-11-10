extern "C"
{
  //! Structure holding shared pointer to an atom.
  struct SetIterator
  { 
    PyObject_HEAD
    Set* ptr_set;
    std::set<std::string>::iterator i_first;
    bool is_first;
  };
  //! Returns Self.
  static PyObject* setiterator_iter(PyObject* _in);
  //! Returns next object.
  static PyObject* setiterator_next(SetIterator* _in);
  //! Function to deallocate a string atom.
  static dealloc setiterator_dealloc(Set *_self);
  //! Creates iterator.
  static PyObject* setiterator_create(Set* _in);
};
// Returns Self.
static PyObject* setiterator_iter(PyObject* _in) { Py_INCREF(_in); return _in; }
// Returns next object.
static PyObject* setiterator_next(SetIterator* _in)
{
  if(_in->i_first != _in->ptr_set->end()) 
  {
    if(_in->is_first) _in->is_first = false; else ++_in->i_first;
    if(_in->i_first != _in->ptr_set->end()) 
      return PyString_FromString(i_first->c_str());
  }
  PyErr_SetNone(PyExc_StopIteration);
  return NULL;
}
//! Function to deallocate a string atom.
static dealloc setiterator_dealloc(SetIterator *_self)
{
  Set* dummy = _self->ptr_set;
  _self->ptr_set = NULL;
  if(dummy) Py_DECREF(dummy);
}

static PyTypeObject setiterator_type = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "atom.SetIterator",        /*tp_name*/
    sizeof(SetIterator),       /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)setiterator_dealloc,   /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,             
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    (reprfunc)set_str,         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_HAVE_ITER,
    "Iterator over atomic occupations of a site.\n\n",
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    (getiterfunc)setiterator_iter,  /* tp_iter */
    (iternextfunc)setiterator_next, /* tp_iternext */
};

// Creates iterator.
static PyObject* setiterator_create(Set* _in)
{
  SetIterator *input = NULL;
  result = PyObject_New(Set, (&setiterator_type));
  if(result == NULL) return NULL;
  result->ptr_set = _in;
  Py_INCREF((PyObject*)result->ptr_set);
  result->i_first = _in->begin();
  result->is_first = true;
  return (PyObject*) result;
}


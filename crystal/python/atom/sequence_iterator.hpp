extern "C"
{
  //! Structure holding shared pointer to an atom.
  struct SequenceIterator
  { 
    PyObject_HEAD
    Sequence* ptr_seq;
    std::vector<std::string>::iterator i_first;
    bool is_first;
  };
  //! Returns Self.
  static PyObject* seqiterator_iter(PyObject* _in);
  //! Returns next object.
  static PyObject* seqiterator_next(SequenceIterator* _in);
  //! Function to deallocate a string atom.
  static void seqiterator_dealloc(SequenceIterator *_self);
  //! Creates iterator.
  static PyObject* sequenceiterator_create(Sequence* _in);
};
// Returns Self.
static PyObject* seqiterator_iter(PyObject* _in) { Py_INCREF(_in); return _in; }
// Returns next object.
static PyObject* seqiterator_next(SequenceIterator* _in)
{
  if(_in->i_first != _in->ptr_seq->ptr_seq->end()) 
  {
    if(_in->is_first) _in->is_first = false; else ++_in->i_first;
    if(_in->i_first != _in->ptr_seq->ptr_seq->end()) 
      return PyString_FromString(_in->i_first->c_str());
  }
  PyErr_SetNone(PyExc_StopIteration);
  return NULL;
}
//! Function to deallocate a string atom.
static void seqiterator_dealloc(SequenceIterator *_self)
{
  Sequence* dummy = _self->ptr_seq;
  _self->ptr_seq = NULL;
  if(dummy) Py_DECREF(dummy);
}

static PyTypeObject sequenceiterator_type = {
    PyObject_HEAD_INIT(NULL)
    0,                                          /*ob_size*/
    "atom.SequenceIterator",                    /*tp_name*/
    sizeof(SequenceIterator),                   /*tp_basicsize*/
    0,                                          /*tp_itemsize*/
    (destructor)seqiterator_dealloc,            /*tp_dealloc*/
    0,                                          /*tp_print*/
    0,                                          /*tp_getattr*/
    0,                                          /*tp_setattr*/
    0,                                          /*tp_compare*/
    0,                                          /*tp_repr*/
    0,                                          
    0,                                          /*tp_as_sequence*/
    0,                                          /*tp_as_mapping*/
    0,                                          /*tp_hash */
    0,                                          /*tp_call*/
    0,                                          /*tp_str*/
    0,                                          /*tp_getattro*/
    0,                                          /*tp_setattro*/
    0,                                          /*tp_as_buffer*/
    Py_TPFLAGS_HAVE_ITER,
    "Iterator over atomic occupations of a site.\n\n",
    0,                                          /* tp_traverse */
    0,                                          /* tp_clear */
    0,                                          /* tp_richcompare */
    0,		                                /* tp_weaklistoffset */
    (getiterfunc)seqiterator_iter,              /* tp_iter */
    (iternextfunc)seqiterator_next,             /* tp_iternext */
};

// Creates iterator.
static PyObject* sequenceiterator_create(Sequence* _in)
{
  SequenceIterator *result = NULL;
  result = PyObject_New(SequenceIterator, (&sequenceiterator_type));
  if(result == NULL) return NULL;
  result->ptr_seq = _in;
  Py_INCREF((PyObject*)result->ptr_seq);
  result->i_first = _in->ptr_seq->begin();
  result->is_first = true;
  return (PyObject*) result;
}


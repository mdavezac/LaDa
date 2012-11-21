struct StructureIterator
{
  PyObject_HEAD;
  PyStructureObject* parent; 
  std::vector<Atom>::const_iterator i_first;
  bool is_first;
};
//! Returns Self.
PyObject* structureiterator_iter(PyObject* _self) { Py_INCREF(_self); return _self; }
//! Returns next object.
PyObject* structureiterator_next(StructureIterator* _self);
//! Function to deallocate a string atom.
void structureiterator_dealloc(StructureIterator *_self);
//! Creates iterator.
PyObject* structureiterator_create(PyStructureObject* _self);


// Returns next object.
PyObject* structureiterator_next(StructureIterator* _self)
{
  if(_self->i_first != _self->parent->atoms.end()) 
  {
    if(_self->is_first) _self->is_first = false; else ++_self->i_first;
    if(_self->i_first != _self->parent->atoms.end()) return _self->i_first->new_ref();
  }
  PyErr_SetNone(PyExc_StopIteration);
  return NULL;
}
//! Function to deallocate a string atom.
void structureiterator_dealloc(StructureIterator *_self)
{
  PyObject* dummy = (PyObject*)_self->parent;
  _self->parent = NULL;
  Py_XDECREF(dummy);
}
// Creates iterator.
PyObject* structureiterator_create(PyStructureObject* _in)
{
  StructureIterator *result = PyObject_New(StructureIterator, structureiterator_type());
  if(result == NULL) return NULL;
  result->parent = _in;
  Py_INCREF(result->parent);
  new(&result->i_first) std::vector<Atom>::iterator(_in->atoms.begin());
  result->is_first = true;
  return (PyObject*) result;
}

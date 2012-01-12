namespace LaDa
{
  namespace crystal
  {
    extern "C" 
    {
      struct StructureIterator
      {
        PyObject_HEAD;
        StructureData* parent; 
        std::vector<Atom>::const_iterator i_first;
        bool is_first;
      };
      //! Returns Self.
      static PyObject* structureiterator_iter(PyObject* _self) { Py_INCREF(_self); return _self; }
      //! Returns next object.
      static PyObject* structureiterator_next(StructureIterator* _self);
      //! Function to deallocate a string atom.
      static void structureiterator_dealloc(StructureIterator *_self);
      //! Creates iterator.
      static PyObject* structureiterator_create(StructureData* _self);
    }

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
    static void structureiterator_dealloc(StructureIterator *_self)
    {
      PyObject* dummy = (PyObject*)_self->parent;
      _self->parent = NULL;
      Py_XDECREF(dummy);
    }
    // Creates iterator.
    static PyObject* structureiterator_create(StructureData* _in)
    {
      StructureIterator *result = PyObject_New(StructureIterator, structureiterator_type());
      if(result == NULL) return NULL;
      result->parent = _in;
      Py_INCREF(result->parent);
      result->i_first = _in->atoms.begin();
      result->is_first = true;
      return (PyObject*) result;
    }
  }
}

namespace LaDa
{
  namespace vff
  {
    extern "C" 
    {
      //! Describes a node in a first neighbor net.
      struct BondIterator
      {
        PyObject_HEAD 
        //! Center of the bonds which are being iterated.
        NodeData* parent;
        //! First and last iterator.
        std::vector<EdgeData*> :: iterator i_first, i_end;
        //! Whether the first next was called. 
        bool is_first; 
        //! Size on creation.
        std::vector<EdgeData*> :: size_type N;
      };
      //! Returns Self.
      static PyObject* bond_iterator_iter(PyObject* _self) { Py_INCREF(_self); return _self; }
      //! Returns next object.
      static PyObject* bond_iterator_next(BondIterator* _self);
      //! Function to deallocate a string atom.
      static void bond_iterator_dealloc(BondIterator *_self);
      //! Creates iterator.
      static PyObject* bond_iterator_create(NodeData* _self);
    }

    // Returns next object.
    PyObject* bond_iterator_next(BondIterator* _self)
    {
      if(_self->N != _self->parent->bonds.size())
      {
        LADA_PYERROR(IndexError, "Bonds changed mid-iteration.");
        return NULL;
      }
      if(_self->i_first != _self->i_end) 
      {
        if(_self->is_first) _self->is_first = false; else ++_self->i_first;
        if(_self->i_first != _self->i_end)
          return edge_to_tuple(_self->parent, *_self->i_first);
      }
      PyErr_SetNone(PyExc_StopIteration);
      return NULL;
    }
    //! Function to deallocate a string atom.
    static void bond_iterator_dealloc(BondIterator *_self)
    {
      PyObject* dummy = (PyObject*)_self->parent;
      _self->parent = NULL;
      Py_XDECREF(dummy);
    }
    //! Creates iterator.
    static PyObject* bond_iterator_create(NodeData* _in)
    {
      BondIterator *result = PyObject_New(BondIterator, bonditerator_type());
      if(result == NULL) return NULL;
      result->parent = _in;
      Py_INCREF(result->parent);
      result->i_first = _in->bonds.begin();
      result->i_end = _in->bonds.end();
      result->is_first = true;
      result->N = _in->bonds.size();
      return (PyObject*) result;
    }
    //! Returns self.
    static PyObject* get_self(PyObject* _self) { Py_INCREF(_self); return _self; }
  }
}

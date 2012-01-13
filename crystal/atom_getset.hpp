namespace LaDa
{
  namespace crystal
  {
    extern "C"
    {
      //! Returns position as a numpy array. 
      static PyObject* lada_atom_getpos(AtomData *_self, void *closure);
      //! Sets position from a sequence of three numbers.
      static int lada_atom_setpos(AtomData *_self, PyObject *_value, void *_closure);
      //! Returns type python object.
      static PyObject* lada_atom_gettype(AtomData *_self, void *closure);
      //! Sets type python object.
      static int lada_atom_settype(AtomData *_self, PyObject *_value, void *_closure);
    }
  
    // Returns position as a numpy array. 
    // Numpy does not implement python's cyclic garbage, hence new wrapper need be
    // created each call.
    static PyObject* lada_atom_getpos(AtomData *_self, void *closure)
    {
      npy_intp dims[1] = {3};
      int const value = math::numpy::type<math::rVector3d::Scalar>::value;
      PyArrayObject* result = (PyArrayObject*) PyArray_SimpleNewFromData(1, dims, value, _self->pos.data());
      if(result == NULL) return NULL;
      result->base = (PyObject*)_self;
      Py_INCREF(_self); // Increfed as base of array.
      return (PyObject*)result;
    }
    // Sets position from a sequence of three numbers.
    static int lada_atom_setpos(AtomData *_self, PyObject *_value, void *_closure)
    {
      if(_value == NULL)
      {
        LADA_PYERROR(TypeError, "Cannot delete pos attribute.");
        return -1;
      }
      return math::convert_to_vector(_value, _self->pos) ? 0: -1;
    }
    // Returns type python object.
    PyObject* lada_atom_gettype(AtomData *_self, void *closure)
      { Py_INCREF(_self->type); return _self->type; }
      
    // Sets type python object.
    int lada_atom_settype(AtomData *_self, PyObject *_value, void *_closure)
    {
      if(_value == NULL)
      {
        LADA_PYERROR(TypeError, "Cannot delete type. You can set it to None though.");
        return -1;
      }
      PyObject* dummy = _self->type;
      _self->type = _value;
      Py_INCREF(_value);
      Py_DECREF(dummy);
      return 0;
    }
  } // namespace Crystal
} // namespace LaDa
  

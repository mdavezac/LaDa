namespace LaDa
{
  namespace crystal
  {
    extern "C"
    {
      //! Returns cell as a numpy array. 
      static PyObject* structure_getcell(StructureData *_self, void *closure);
      //! Sets cell from a sequence of 3x3 numbers.
      static int structure_setcell(StructureData *_self, PyObject *_value, void *_closure);
      // Returns the scale.
      static PyObject* structure_getscale(StructureData *_self, void *closure);
      //! Sets the scale from a number.
      static int structure_setscale(StructureData *_self, PyObject *_value, void *_closure);
    }
  
    // Returns cell as a numpy array. 
    // Numpy does not implement python's cyclic garbage, hence new wrapper need be
    // created each call.
    static PyObject* structure_getcell(StructureData *_self, void *closure)
    {
      npy_intp dims[2] = {3, 3};
      int const value = math::numpy::type<math::rMatrix3d::Scalar>::value;
      PyArrayObject *result = (PyArrayObject*) PyArray_SimpleNewFromData(2, dims, value, _self->cell.data());
      if(result == NULL) return NULL;
      result->base = (PyObject*)_self;
      Py_INCREF(_self); // Increfed as base of array.
      return (PyObject*)result;
    }
    // Sets cell from a sequence of three numbers.
    static int structure_setcell(StructureData *_self, PyObject *_value, void *_closure)
    {
      if(_value == NULL)
      {
        LADA_PYERROR(TypeError, "Cannot delete cell attribute.");
        return -1;
      }
      return math::convert_to_cell(_value, _self->cell) ? 0: -1;
    }
  
    // Returns the scale of the structure.
    static PyObject* structure_getscale(StructureData *_self, void *closure)
      { return PyFloat_FromDouble(_self->scale); }
    // Sets the scale of the structure from a number.
    static int structure_setscale(StructureData *_self, PyObject *_value, void *_closure)
    {
      if(_value == NULL) 
      {
        LADA_PYERROR(TypeError, "Cannot delete scale attribute.");
        return -1;
      }
      if(PyFloat_Check(_value)) _self->scale = PyFloat_AS_DOUBLE(_value); 
      else if(PyInt_Check(_value)) _self->scale = PyInt_AS_LONG(_value); 
      else if(python::Quantity::isinstance(_value))
        try { _self->scale = python::Quantity(_value).get("angstrom"); }
        catch(...) { return -1; }
      else
      {
        LADA_PYERROR(TypeError, "Input to scale is not an acceptable type.");
        return -1;
      }
      return 0;
    }
  }
} 

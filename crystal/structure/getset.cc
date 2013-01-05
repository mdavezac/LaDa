//! Returns cell as a numpy array. 
PyObject* structure_getcell(PyStructureObject *_self, void *closure)
   { return python::numpy::wrap_to_numpy(_self->cell, (PyObject*)_self); }
//! Sets cell from a sequence of 3x3 numbers.
int structure_setcell(PyStructureObject *_self, PyObject *_value, void *_closure);
// Returns the scale.
PyObject* structure_getscale(PyStructureObject *_self, void *closure);
//! Sets the scale from a number.
int structure_setscale(PyStructureObject *_self, PyObject *_value, void *_closure);
//! Gets the volume of the structure
PyObject* structure_getvolume(PyStructureObject *_self, void *_closure)
{
  types::t_real const scale = python::get_quantity(_self->scale);
  types::t_real const result = std::abs(_self->cell.determinant() * std::pow(scale, 3));
  return python::fromC_quantity(result, _self->scale);
}


// Sets cell from a sequence of three numbers.
int structure_setcell(PyStructureObject *_self, PyObject *_value, void *_closure)
{
  if(_value == NULL)
  {
    PYLADA_PYERROR(TypeError, "Cannot delete cell attribute.");
    return -1;
  }
  return python::numpy::convert_to_matrix(_value, _self->cell) ? 0: -1;
}

// Returns the scale of the structure.
PyObject* structure_getscale(PyStructureObject *_self, void *closure)
{
  Py_INCREF(_self->scale); 
  return _self->scale;
}
// Sets the scale of the structure from a number.
int structure_setscale(PyStructureObject *_self, PyObject *_value, void *_closure)
{
  if(_value == NULL) 
  {
    PYLADA_PYERROR(TypeError, "Cannot delete scale attribute.");
    return -1;
  }
  PyObject *result = python::fromPy_quantity(_value, _self->scale);
  if(not result) return -1;
  PyObject *dummy = _self->scale;
  _self->scale = result;
  Py_DECREF(dummy);
  return 0;
}

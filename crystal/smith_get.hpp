namespace LaDa
{
  namespace crystal
  {
    extern "C"
    {
      //! Returns transformation matrix as numpy array.
      static PyObject* smithtransform_getmatrix(SmithTransformData *_self, void *closure)
         { return python::wrap_to_numpy(_self->matrix, (PyObject*)_self); }
      //! Returns periodicity vector as numpy array.
      static PyObject* smithtransform_getvector(SmithTransformData *_self, void *closure)
         { return python::wrap_to_numpy(_self->vector, (PyObject*)_self); }
      //! Number of unit-cells in the supercell.
      static PyObject* smithtransform_getsize(SmithTransformData *_self, void *closure)
         { return python::wrap_to_numpy(_self->vector, (PyObject*)_self); }
    }
  }
} 

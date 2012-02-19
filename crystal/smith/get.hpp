namespace LaDa
{
  namespace crystal
  {
    extern "C"
    {
      //! Returns transformation matrix as numpy array.
      static PyObject* smithtransform_gettransform(SmithTransformData *_self, void *closure)
         { return python::wrap_to_numpy(_self->transform, (PyObject*)_self); }
      //! Returns periodicity quotient as numpy array.
      static PyObject* smithtransform_getquotient(SmithTransformData *_self, void *closure)
         { return python::wrap_to_numpy(_self->quotient, (PyObject*)_self); }
      //! Number of unit-cells in the supercell.
      static PyObject* smithtransform_getsize(SmithTransformData *_self, void *closure)
         { return python::wrap_to_numpy(_self->quotient, (PyObject*)_self); }
    }
  }
} 

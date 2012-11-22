//! Returns transformation matrix as numpy array.
PyObject* hftransform_gettransform(PyHFTObject *_self, void *closure)
   { return python::numpy::wrap_to_numpy(_self->transform, (PyObject*)_self); }
//! Returns periodicity quotient as numpy array.
PyObject* hftransform_getquotient(PyHFTObject *_self, void *closure)
   { return python::numpy::wrap_to_numpy(_self->quotient, (PyObject*)_self); }
//! Number of unit-cells in the supercell.
PyObject* hftransform_getsize(PyHFTObject *_self, void *closure)
   { return PyLong_FromLong(   _self->quotient[0] 
                             * _self->quotient[1]
                             * _self->quotient[2] ); }

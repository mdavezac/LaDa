
namespace LaDa
{
  namespace crystal
  {
    extern "C" 
    {
      //! Returns a deepcopy of the atom.
      static PyObject* smithtransform_copy(SmithTransformData* _self)
        { return (PyObject*) PySmithTransform_Copy(_self, NULL); }
      //! Implements deepcopy.
      static PyObject* smithtransform_deepcopy(SmithTransformData* _self, PyObject* _memo)
        { return (PyObject*) PySmithTransform_Copy(_self, _memo); }
      //! Implements shallow copy.
      static PyObject* smithtransform_shallowcopy(SmithTransformData* _self)
        { Py_INCREF(_self); return (PyObject*)_self; }
      //! Implements getstate for pickling.
      static PyObject* smithtransform_getstate(SmithTransformData* _self);
      //! Implements setstate for pickling.
      static PyObject* smithtransform_setstate(SmithTransformData* _self, PyObject *_dict);
      //! Implements reduce for pickling.
      static PyObject* smithtransform_reduce(SmithTransformData* _self);
      //! Computes Z-group indices of position \a _pos.
      static PyObject* smithtransform_indices(SmithTransformData* _self, PyObject* _args);
      // Computes flat smith index from non-flat smith index.
      static PyObject* smithtransform_flatten_indices( SmithTransformData* _self,
                                                       PyObject* _args, PyObject *_kwargs );
      //! Computes flat index into Z-group from atomic position.
      static PyObject* smithtransform_flat_index( SmithTransformData* _self,
                                                  PyObject* _args, PyObject *_kwargs );
    }

    // Implements __reduce__ for pickling.
    PyObject* smithtransform_reduce(SmithTransformData* _self)
    {
      // Creates return tuple of three elements.
      python::Object type = PyObject_Type((PyObject*)_self);
      if(not type) return NULL;
      // Second element is a null tuple, argument to the callable type above.
      python::Object tuple = PyTuple_New(0);
      if(not tuple) return NULL;
      // Third element is the state of this object.
      char getstate[] = "__getstate__";
      python::Object state = PyObject_CallMethod((PyObject*)_self, getstate, NULL, NULL);
      if(not state) return NULL;

      return PyTuple_Pack(3, type.borrowed(), tuple.borrowed(), state.borrowed());
    }

    // Implements getstate for pickling.
    PyObject* smithtransform_getstate(SmithTransformData* _self)
    {
      // get cell attribute.
      python::Object cell = smithtransform_gettransform(_self, NULL);
      if(not cell) return NULL;
      // get scale attribute.
      python::Object scale = smithtransform_getquotient(_self, NULL);
      if(not scale) return NULL;
      return PyTuple_Pack(2, cell.borrowed(), scale.borrowed());
    }

    // Implements setstate for pickling.
    PyObject* smithtransform_setstate(SmithTransformData* _self, PyObject *_tuple)
    {
      if(not PyTuple_Check(_tuple))
      {
        LADA_PYERROR(TypeError, "Expected state to be a tuple.");
        return NULL;
      }
      if(PyTuple_Size(_tuple) != 2)
      {
        LADA_PYERROR(TypeError, "Expected state to be a 2-tuple.");
        return NULL;
      }
      // first cell and scale.
      if(not python::convert_to_matrix(PyTuple_GET_ITEM(_tuple, 0), _self->transform)) return NULL;
      if(not python::convert_to_vector(PyTuple_GET_ITEM(_tuple, 1), _self->quotient)) return NULL;
      Py_RETURN_NONE;
    }

    // defines macros also used in smith.
#   include "smith_macro.hpp"
    // Computes flat smith index from non-flat smith index.
    static PyObject* smithtransform_flatten_indices( SmithTransformData* _self,
                                                     PyObject* _args, PyObject *_kwargs )
    {
      PyObject *posatom = NULL;
      int site = -1;
      static char *kwlist[] = { const_cast<char*>("pos"), const_cast<char*>("site"), NULL };
      if(not PyArg_ParseTupleAndKeywords(_args, _kwargs, "O|i:pos_index", kwlist, &posatom, site) )
        return NULL;
      math::iVector3d pos;
      if(not python::convert_to_vector(posatom, pos)) return NULL;
      LADA_SMITHTRANSFORM_SHARED0(_self->quotient, pos, site);
      return PyInt_FromLong(flat_result);
    }
    //! Computes flat smith index from non-flat smith index, including sites.
    static PyObject* smithtransform_flat_index( SmithTransformData* _self,
                                                PyObject* _args, PyObject *_kwargs )
    {
      PyObject *posatom = NULL;
      int site = -1;
      static char *kwlist[] = { const_cast<char*>("indices"), const_cast<char*>("site"), NULL };
      if(not PyArg_ParseTupleAndKeywords(_args, _kwargs, "O|i:flat_index", kwlist, &posatom, &site) )
        return NULL;
      math::rVector3d pos;
      if(not python::convert_to_vector(posatom, pos)) return NULL;
      LADA_SMITHTRANSFORM_SHARED1(_self->quotient, _self->transform, pos, LADA_PYERROR, return NULL);
      LADA_SMITHTRANSFORM_SHARED0(_self->quotient, vector_result, site);
      return PyInt_FromLong(flat_result);
    }
    // Computes smith indices of position \a _pos.
    static PyObject* smithtransform_indices(SmithTransformData* _self, PyObject* _args)
    {
      math::rVector3d pos;
      if(not python::convert_to_vector(_args, pos)) return NULL;
      LADA_SMITHTRANSFORM_SHARED1(_self->quotient, _self->transform, pos, LADA_PYERROR, return NULL);
      return python::wrap_to_numpy(vector_result);
    }
    // undefs macros also used in smith.
#   include "smith_macro.hpp"
  }
}


extern "C" 
{
  //! Returns size of structure.
  static Py_ssize_t LADA_NAME(size)(LADA_TYPE* _self) { return _self->structure->atoms.size(); }
  //! Returns size of structure.
  static PyObject* LADA_NAME(getitem)(LADA_TYPE* _self, PyObject *_key);
}

static PyObject* LADA_NAME(getitem)(LADA_TYPE* _self, PyObject *_key)
{
  if(PyInt_Check(_key)) 
  {
    int i = PyInt_AsLong(_key);
    if(i < 0) i += _self->structure->atoms.size();
    if(i < 0 or i >= _self->structure->atoms.size()) 
    {
      LADA_PYERROR(IndexError, "Index out of range.");
      return NULL;
    }
    return _self->structure->atoms[i].pyself();
  }
  std::cout << _key->ob_type->tp_name << "\n";
  LADA_PYERROR(TypeError, "Could not make sense of index to atom.");
  return NULL;
}

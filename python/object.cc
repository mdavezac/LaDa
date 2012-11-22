void object_reset(PyObject*& _object, PyObject *_in)
{
  PyObject * const dummy(_object);
  _object = _in;
  Py_XINCREF(_object);
  Py_XDECREF(dummy);
}

bool object_equality_op(Object const& _self, Object const &_b)
{
  if(not _self.hasattr("__eq__"))
    BOOST_THROW_EXCEPTION( error::TypeError() << 
                           error::string("No __eq__ member function found "
                                         "when comparing objects.") );
  Object methodname = PyString_FromString("__eq__");
  if(not methodname)
    BOOST_THROW_EXCEPTION(error::internal() << 
                          error::string("Could not create string."));
  Object result = PyObject_CallMethodObjArgs(_self.borrowed(), 
                                             methodname.borrowed(), _b.borrowed(), NULL);
  if(not result) 
    BOOST_THROW_EXCEPTION(error::internal() <<
                          error::string("Python exception thrown when comparing objects."));
  // Try reflected operation.
  if(result.borrowed() == Py_NotImplemented)
  {
    if(not _b.hasattr("__eq__"))
    {
      if(_b.borrowed()->ob_type != _self.borrowed()->ob_type) return false;
      BOOST_THROW_EXCEPTION( error::TypeError() << 
                             error::string( "No implementation of "
                                            "equality between these "
                                            "two object has been found."));
    }
    result.reset( PyObject_CallMethodObjArgs( _b.borrowed(), 
                                              methodname.borrowed(), 
                                              _self.borrowed(), NULL) );
    if(not result)
      BOOST_THROW_EXCEPTION(error::TypeError() << 
                            error::string("Python exception thrown when comparing objects."));
    if(result.borrowed() == Py_NotImplemented)
    {
      if(_b.borrowed()->ob_type != _self.borrowed()->ob_type) return false;
      BOOST_THROW_EXCEPTION( error::TypeError() <<
                             error::string( "No implementation of equality between"
                                            " these two object has been found.") );
    }
  }
  if(PyBool_Check(result.borrowed())) return result.borrowed() == Py_True;
  if(PyInt_Check(result.borrowed())) return PyInt_AS_LONG(result.borrowed()) != 0;
  BOOST_THROW_EXCEPTION( error::ValueError() <<
                         error::string("Could not make sense of return "
                                       "of comparison function."));
};

std::ostream& operator<< (std::ostream &stream, Object const &_ob)
{
  PyObject* const repr = PyObject_Repr(_ob.borrowed());
  if(not repr) BOOST_THROW_EXCEPTION(error::internal());
  char const * const result = PyString_AS_STRING(repr);
  if(not result) BOOST_THROW_EXCEPTION(error::internal()); 
  return stream << result;
}

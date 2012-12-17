void object_reset(PyObject*& _object, PyObject *_in)
{
  PyObject * const dummy(_object);
  _object = _in;
  Py_XINCREF(_object);
  Py_XDECREF(dummy);
}

bool object_equality_op(Object const& _self, Object const &_b)
{
  if(not _self)
    BOOST_THROW_EXCEPTION( error::internal()
                           << error::string("Object yet uninitialized.") );
  if(not _b)
    BOOST_THROW_EXCEPTION( error::internal()
                           << error::string("Object yet uninitialized.") );

  PyObject *globals = PyEval_GetBuiltins();
  if(not globals)
    BOOST_THROW_EXCEPTION( error::internal() 
                           << error::string("Could not get builtins.") );
  python::Object locals = PyDict_New();
  if(not locals)
    BOOST_THROW_EXCEPTION( error::internal() 
                           << error::string("Could not create local dict.") );
  if(PyDict_SetItemString(locals.borrowed(), "a", _self.borrowed()) < 0)
    BOOST_THROW_EXCEPTION( error::internal() 
                           << error::string("Could not set item in local dict.") );
  if(PyDict_SetItemString(locals.borrowed(), "b", _b.borrowed()) < 0)
    BOOST_THROW_EXCEPTION( error::internal() 
                           << error::string("Could not set item in local dict.") );
  
  python::Object const result(PyRun_String( "(a == b) == True", Py_eval_input, globals, 
                                            locals.borrowed() ));
  return result.borrowed() == Py_True;
};

std::ostream& operator<< (std::ostream &stream, Object const &_ob)
{
  PyObject* const repr = PyObject_Repr(_ob.borrowed());
  if(not repr) BOOST_THROW_EXCEPTION(error::internal());
  char const * const result = PyString_AS_STRING(repr);
  if(not result) BOOST_THROW_EXCEPTION(error::internal()); 
  return stream << result;
}

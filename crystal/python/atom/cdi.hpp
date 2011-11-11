extern "C" 
{ 
  //! Function to deallocate a string atom.
  static void atomstr_Dealloc(AtomStr *_self);
  //! Function to allocate a string atom.
  static PyObject* atomstr_new(PyTypeObject *_type, PyObject *_args, PyObject *_kwargs);
  //! Function to initialize a string atom.
  static int atomstr_init(AtomStr* _self, PyObject* _args, PyObject *_kwargs);
  //! Traverses to back-reference.
  static int traverse(AtomStr *_self, visitproc _visit, void *_arg);
  //! Clears back reference.
  static int gcclear(AtomStr *_self);
}

// Function to deallocate a string atom.
static void atomstr_Dealloc(AtomStr *_self)
{
  if(_self->weakreflist != NULL)
    PyObject_ClearWeakRefs((PyObject *) _self);

  gcclear(_self);
  boost::shared_ptr< LaDa::crystal::AtomData< std::string > > dummy;
  dummy.swap(_self->atom);
  _self->ob_type->tp_free((PyObject*)_self);
}

//! Function to allocate a string atom.
static PyObject* atomstr_new(PyTypeObject *_type, PyObject *_args, PyObject *_kwargs)
{
  PyObject *const pydict = PyDict_New();
  if(pydict == NULL) return NULL;
  AtomStr *self;
  self = (AtomStr*)_type->tp_alloc(_type, 0);
  if(self == NULL) return NULL;
  self->weakreflist = NULL;
  boost::shared_ptr< crystal::AtomData<std::string> > dummy(new LaDa::crystal::AtomData<std::string>);
  if(not dummy)
  {
    Py_DECREF(self);
    PyErr_SetString( PyException<error::internal>::exception().ptr(), 
                     "Could not create atom.\n" );
    return NULL;
  } 
  self->atom.swap(dummy);
  self->dictionary = pydict;

  // create numpy array referencing position.
  npy_intp dims[1] = {3};
  int const value = math::numpy::type<math::rVector3d::Scalar>::value;
  self->position = (PyArrayObject*) PyArray_SimpleNewFromData(1, dims, value, self->atom->pos.data());
  bool const err = PyErr_Occurred() != NULL;
  if(self->position == NULL or err)
  {
    Py_DECREF(self);
    if(self->position != NULL) Py_DECREF(self->position);
    if(err)
      PyErr_SetString( PyException<error::internal>::exception().ptr(), 
                       "Could not create position array in atom.\n" );
    return NULL;
  } 
  self->position->base = (PyObject*)self;
  Py_INCREF(self);

  return (PyObject*) self;
}

// Function to initialize a string atom.
static int atomstr_init(AtomStr* _self, PyObject* _args, PyObject *_kwargs)
{
  try
  {
    bp::dict kwargs;
    if(_kwargs) kwargs = bp::dict(bp::handle<>(bp::borrowed(_kwargs)));
    bp::object args(bp::handle<>(bp::borrowed(_args)));
    if(bp::len(args) == 0 and bp::len(kwargs) == 0) return 0;
    
    // first looks to identify position in input argument or dictionary..
    bool found_position = false;
    if( bp::len(args) >= 1 and lp::is_position(args[0]) ) 
    {
      lp::extract_position(args[0], _self->atom->pos);
      found_position = true;
      args = args.slice(1, bp::slice_nil());
    }
    else if( bp::len(args) >= 3)
    {
      if(not lp::is_position(args))
      {
        PyErr_SetString( PyException<error::TypeError>::exception().ptr(),
                         "First three arguments could not be translated to a position." );
        return -1;
      }
      lp::extract_position(args.slice(0, 3), _self->atom->pos);
      found_position = true;
      args = args.slice(3, bp::slice_nil());
    }
    if( kwargs.has_key("position") )
    {
      if(found_position)
      {
        PyErr_SetString( PyException<error::TypeError>::exception().ptr(),
                         "Multiple value for position." );
        return -1;
      }
      lp::extract_position(kwargs["position"], _self->atom->pos);
      PyDict_DelItemString(kwargs.ptr(), "position");
      found_position = true;
    }
    
    // Now looks for specie.
    bool found_specie = false;
    if( bp::len(args) > 1) 
    {
      PyErr_SetString( PyException<error::TypeError>::exception().ptr(),
                       "Did not understand argument. Did you input more than one specie?" );
      return -1;
    }
    else if( bp::len(args) == 1 )
    {
      if(not lp::is_specie<std::string>(args[0]))
      {
        PyErr_SetString( PyException<error::TypeError>::exception().ptr(),
                         "Argument did not translate to a type." );
        return -1;
      }
      lp::extract_specie(args[0], _self->atom->type);
      found_specie = true;
    }
    if( kwargs.has_key("type") )
    {
      if(found_specie)
      {
        PyErr_SetString( PyException<error::TypeError>::exception().ptr(),
                         "Multiple value for specie." );
        return -1;
      }
      lp::extract_specie(kwargs["type"], _self->atom->type);
      PyDict_DelItemString(kwargs.ptr(), "type");
    }
    // Now freeze and site
    if(kwargs.has_key("site"))
    {
      _self->atom->site = bp::extract<types::t_int>(kwargs["site"]);
      PyDict_DelItemString(kwargs.ptr(), "site");
    }
    if(kwargs.has_key("freeze"))
    {
      _self->atom->freeze = bp::extract<types::t_int>(kwargs["freeze"]);
      PyDict_DelItemString(kwargs.ptr(), "freeze");
    }
    // Now additional attributes.
    PyObject *__dict__ = PyObject_GetAttrString((PyObject*) _self, "__dict__");
    if(__dict__ == NULL)
    {
      PyErr_SetString( PyException<error::internal>::exception().ptr(),
                       "Could not extract __dict__ in atom." );
      return -1;
    }
    int const result = PyDict_Merge(__dict__, kwargs.ptr(), 1);
    Py_DECREF(__dict__); 
    return result;
  }
  catch(bp::error_already_set &e) { return -1; }
}

static int traverse(AtomStr *self, visitproc visit, void *arg)
{
  Py_VISIT(self->position);
  Py_VISIT(self->position->base);
  return 0;
}

static int gcclear(AtomStr *_self)
{ 
  PyArrayObject* position = _self->position;
  _self->position = NULL;
  Py_XDECREF(position);
  return 0;
}

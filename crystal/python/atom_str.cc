#include "LaDaConfig.h"

#include <Python.h>
#include <boost/python/tuple.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/borrowed.hpp>

#include <math/python/python.hpp>
#include <python/exceptions.h>

#include "../atom.h"

using namespace bp = boost::python;
using namespace lp = LaDa::python;

//! Structure holding shared pointer to an atom.
extern "C" struct AtomStr
{
  PyObject_HEAD
  boost::shared_ptr< LaDa::crystal::AtomData< std::string > > atom;
};

//! Function to deallocate a string atom.
extern "C" static AtomStr_Dealloc(AtomStr *_self)
{
  boost::shared_ptr< LaDa::crystal::AtomData< std::string > > dummy;
  dummy.swap(_self->atom);
  _self->ob_type->tp_free((PyObject*)_self);
}

//! Function to allocate a string atom.
extern "C" static PyObject* AtomStr_new(PyTypeObject *_type, PyObject *_args, PyObject *_kwargs)
{
  AtomStr *self;
  self = (AtomStr*)_type->tp_alloc(type, 0);
  self->atom.reset(new LaDa::crystal::AtomData<std::string>);
  if(not self->atom)
  {
    Py_DECREF(self);
    return NULL;
  } 
  self->atom->set_self(self);
  return (PyObject*) self;
}

//! Function to initialize a string atom.
extern "C" static int AtomStr_init(AtomStr* _self, PyObject* _args, PyObject *_kwargs)
{
  bp::dict kwargs(bp::handle<>(bp::borrowed(_kwargs)));
  bp::tuple args(bp::handle<>(bp::borrowed(_args)));

  // first looks to identify position in input argument or dictionary..
  bool found_position = false
  if( bp::len(args) >= 1 and lp::is_position(args[0]) ) 
  {
    lp::extract_position(args[0], _self->atom->pos);
    found_position = true;
    args = args.slice(1, bp::slice_nil);
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
    args = args.slice(3, bp::slice_nil);
  }
  if( kwargs.has_key("position") )
  {
    if(found_position)
    {
      PyErr_SetString( PyException<error::TypeError>::exception().ptr(),
                       "Multiple value for position." )
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
    if(not lp::is_specie(args[0]))
    {
      PyErr_SetString( PyException<error::TypeError>::exception().ptr(),
                       "Argument did not translate to a type." )
      return -1;
    }
    lp::extract(args[0], _self->atom->type);
    found_specie = true;
  }
  if( kwargs.has_dict("type") )
  {
    if(found_specie)
    {
      PyErr_SetString( PyException<error::TypeError>::exception().ptr(),
                       "Multiple value for specie." )
      return -1;
    }
    lp::extract_position(kwargs["type"], _self->atom->type);
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
    _self->atom->freeze = bp::extract<types::t_int>(dict["freeze"]);
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
  int const result = PyDict_Merg(__dict__, kwargs.ptr(), 1);
  Py_DECREF(__dict__); 
  return result;
}

extern "C" static PyObject* AtomStr_getpos(AtomStr *_self, void *closure)
{
  npy_intp dims[1] = {n};
  PyObject *result = PyArray_SimpleNewFromData(1, dims, type<T>::value, _ptr);
  if( result == NULL or PyErr_Occurred() != NULL ) return NULL;
  return result;
}
extern "C" static int AtomStr_setpos(AtomStr *_self, PyObject *_value, void *_closure)
{
  bp::object pos(bp::handle<>(bp::borrowed(_value)));
  if(not is_position(_pos)) 
  {
    PyException<error::TypeError>::exception().ptr()("Input could not be converted to position.");
    return -1;
  }
  extract_position(pos, self->atom->pos);
  return 0;
}

extern "C" static PyObject* AtomStr_getsite(AtomStr *_self, void *closure)
  { return PyInt_FromLong(_self->atom->site); }
extern "C" static int AtomStr_setsite(AtomStr *_self, PyObject *_value, void *_closure)
{
  long const result = PyInt_AsLong(_value);
  if(result == -1 and PyErr_Occurred() != NULL) return -1;
  _self->atom->site = result;
  return 0;
}
extern "C" static PyObject* AtomStr_getfreeze(AtomStr *_self, void *closure)
{
  long result = _self->atom->freeze;
  return PyInt_FromLong(result);
}
extern "C" static int AtomStr_setfreeze(AtomStr *_self, PyObject *_value, void *_closure)
{
  long const result = PyInt_AsLong(_value);
  if(result == -1 and PyErr_Occurred() != NULL) return -1;
  if(result < 0)
  {
    PyErr_SetString( PyException<error::ValueError>::exception().ptr(),  
                     "Cannot set freeze to a negative value." );
    return -1;
  }
  _self->atom->site = result;
  return 0;
}

static PyMethodDef AtomStr_getsetters[] = {
    {"pos", (getter)AtomStr_getpos, (setter)AtomStr_setpos,
     "Position in cartesian coordinates.", NULL},
    {"site", (getter)AtomStr_getsite, (setter)AtomStr_setsite,
     "Site index (integer).", NULL},
    {"freeze", (getter)AtomStr_getfreeze, (setter)AtomStr_setfreeze,
     "Mask to freeze position or type (unsigned integer).", NULL},
    {NULL}  /* Sentinel */
};

static PyTypeObject AtomStrType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "atom.AtomStr",            /*tp_name*/
    sizeof(AtomStr),             /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)AtomStr_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "Atom for which the type is specified as a strings.\n\n"
      "__init__ accepts different kind of input.\n"
      "  - The position can be given as:\n" 
      "      - the first positional argument, in which case "
               "it should be a sequence of three floating points\n"
      "      - the first *three* positional argument\n"
      "      - as a keyword argument ``position``, \n"
      "      - not at all, in which case it default to the origin.\n"
      "  - The type can be given a:\n"
      "      - the first argument if the position is not given or given a keyword.\n"
      "      - the first (and last) argument following the position "
               "if the position is not given as a keyword.\n"
      "      - as a keyword argument ``type``in which case "
               "it should be a string\n"
      "      - not at all, in which case it default to the an empty string.\n"
      "  - The site index and the ``freeze`` parameter can only be given as keywords.\n"
      "  - All other keyword arguments become attributes. "
           "In other words, one could add ``magnetic=0.5`` if one wanted to "
           "specify the magnetic moment of an atom.\n",
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    0,                         /* tp_methods */
    0,                         /* tp_members */
    AtomStr_getsetters,        /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)AtomStr_init,      /* tp_init */
    0,                         /* tp_alloc */
    AtomStr_new,                 /* tp_new */
};

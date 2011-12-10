static PyMethodDef LADA_NAME(methods)[] = {
    {"copy", (PyCFunction)LADA_NAME(copy), METH_NOARGS, "Returns a deepcopy of the atom." },
    {"cast", (PyCFunction)LADA_NAME(cast), METH_VARARGS,
             "If a string atom, returns a sequence atom, and vice-versa.\n\n"
             ":Parameters:\n"
             "  sep : string\n"
             "   Separator between atomic species. Defaults to a comma. "
             "When converting to a sequence from a string, this separator "
             "will split the atomic species. It is a regex in the boost::regex format. "
             "When going from sequences to strings, the separator is placed between species."  },
    {"to_dict", (PyCFunction)LADA_NAME(to_dict), METH_NOARGS,
                "Returns a dictionary with shallow copies of items." },
    {"__copy__", (PyCFunction)LADA_NAME(shallowcopy), METH_NOARGS, "Shallow copy of an atom." },
    {"__deepcopy__", (PyCFunction)LADA_NAME(deepcopy), METH_O, "Deep copy of an atom." },
    {"__getstate__", (PyCFunction)LADA_NAME(getstate), METH_NOARGS, "Implements pickle protocol." },
    {"__setstate__", (PyCFunction)LADA_NAME(setstate), METH_O, "Implements pickle protocol." },
    {"__reduce__", (PyCFunction)LADA_NAME(reduce), METH_NOARGS, "Implements pickle protocol." },
    {NULL}  /* Sentinel */
};

PyTypeObject LADA_NAME(type) = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
#   if LADA_ATOM_NUMBER == 0
      "lada.crystal.cppwrappers.AtomStr",      /*tp_name*/
#   elif LADA_ATOM_NUMBER == 1
      "lada.crystal.cppwrappers.AtomSequence", /*tp_name*/
#   endif
    sizeof(LADA_TYPE),         /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)LADA_NAME(dealloc), /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    (reprfunc)LADA_NAME(repr), /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    LADA_NAME(getattro),       /*tp_getattro*/
    LADA_NAME(setattro),       /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC, /*tp_flags*/
#   if LADA_ATOM_NUMBER == 0
      "Atom for which the type is specified as a string.\n\n"
#   elif LADA_ATOM_NUMBER == 1
      "Atom for which the type is specified as a list of strings.\n\n"
#   endif
        "__init__ accepts different kind of input.\n"
        "  - The position can be given as:\n" 
        "      - the first positional argument, in which case "
                 "it should be a sequence of three floating points\n"
        "      - the first *three* positional argument\n"
        "      - as a keyword argument ``position``, \n"
        "      - not at all, in which case it default to the origin.\n"
#       if LADA_ATOM_NUMBER == 0
          "  - The type can be given as:\n"
          "      - the first argument if the position is not given or given a keyword.\n"
          "      - the first (and last) argument following the position "
                   "if the position is not given as a keyword.\n"
          "      - as a keyword argument ``type``in which case "
                   "it should be a string.\n"
          "      - not at all, in which case it default to the an empty string.\n"
#       elif LADA_ATOM_NUMBER == 1
          "  - The type can be given as:\n"
          "      - the first argument if the position is not given or given a keyword. "
                   "In that case, it should be either a single string or a list of strings.\n"
          "      - the first few arguments if the position is not given or given a keyword. "
                   "In that case, each should be a string.\n"
          "      - the only argument following the position "
                   "In that case it should either a single string or a list of strings.\n"
          "      - A number of argument following the position. "
                   "In that case each argument should be a single string.\n"
          "      - as a keyword argument ``type``in which case "
                   "it should be a string or a list of strings.\n"
          "      - not at all, in which case it default to the an empty list of string.\n"
#       endif
        "  - The site index and the ``freeze`` parameter can only be given as keywords.\n"
        "  - All other keyword arguments become attributes. "
             "In other words, one could add ``magnetic=0.5`` if one wanted to "
             "specify the magnetic moment of an atom.\n",
    (traverseproc)LADA_NAME(traverse),    /* tp_traverse */
    (inquiry)LADA_NAME(gcclear),          /* tp_clear */
    0,		               /* tp_richcompare */
    offsetof(LADA_TYPE, weakreflist), /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    LADA_NAME(methods),        /* tp_methods */
    0,                         /* tp_members */
    LADA_NAME(getsetters),        /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset, dictionary is not reached this way. Must be null. */
    (initproc)LADA_NAME(init),    /* tp_init */
    0,                            /* tp_alloc */
    LADA_NAME(new),               /* tp_new */
};

static PyMappingMethods LADA_NAME(as_mapping)[] = {
   (lenfunc)LADA_NAME(size), 
   (binaryfunc) LADA_NAME(getitem), 
   NULL 
};

static PyMethodDef LADA_NAME(methods)[] = {
    {"copy", (PyCFunction)LADA_NAME(copy), METH_NOARGS, "Returns a deepcopy of the structure." },
    {"cast", (PyCFunction)LADA_NAME(cast), METH_VARARGS,
             "If a string structure, returns a sequence structure, and vice-versa.\n\n"
             ":Parameters:\n"
             "  sep : string\n"
             "   Separator between atomic species. Defaults to a comma." },
    {"to_dict", (PyCFunction)LADA_NAME(to_dict), METH_NOARGS,
                "Returns a dictionary with shallow copies of items." },
    {"__copy__", (PyCFunction)LADA_NAME(shallowcopy), METH_NOARGS, "Shallow copy of an atom." },
    {"__deepcopy__", (PyCFunction)LADA_NAME(deepcopy), METH_O, "Deep copy of an atom." },
    {"__getstate__", (PyCFunction)LADA_NAME(getstate), METH_NOARGS, "Implements pickle protocol." },
    {"__setstate__", (PyCFunction)LADA_NAME(setstate), METH_O, "Implements pickle protocol." },
    {"__reduce__", (PyCFunction)LADA_NAME(reduce), METH_NOARGS, "Implements pickle protocol." },
    {"add_atom", (PyCFunction)LADA_NAME(add_atom), METH_KEYWORDS,
                 "Adds an atom to a structure.\n"
                 "The parameters can be either an atom, a reference (not a copy) too which "
                 "is added to the structure, or it accepts the same arguments and keywords "
                 "as `Atom.__init__`" },
    {NULL}  /* Sentinel */
};

PyTypeObject LADA_NAME(type) = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
#   if LADA_ATOM_NUMBER == 0
      "lada.crystal.cppwrappers.StructureStr",      /*tp_name*/
#   elif LADA_ATOM_NUMBER == 1
      "lada.crystal.cppwrappers.StructureSequence", /*tp_name*/
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
    LADA_NAME(as_mapping),     /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    LADA_NAME(getattro),       /*tp_getattro*/
    LADA_NAME(setattro),       /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC, /*tp_flags*/
#   if LADA_ATOM_NUMBER == 0
      "Structure for which the type is specified as a string.\n\n"
#   elif LADA_ATOM_NUMBER == 1
      "Structure for which the type is specified as a list of strings.\n\n"
#   endif
      ":Parameters:\n"
      "  cell : 3x3 sequence of floats\n"
      "    Cell vectors in units of ``scale`` and cartesian coordinates. "
           "It should be given in the form of a matrix. eg ``cell[:, 0] = a0``, "
           "with ``a0`` the first cell vector. Note that most First Principles "
           "code require the cell vectors as rows rather than columns, as expected here. "
           "Blame the use of Fortran versus more natural languages."
           "It is expected, though not enforced, that ``structure.scale * structure.cell`` "
           "is in angstrom and cartesian coordinates. Eventually, such an arbitrary decision "
           "should lead to fewer mistakes.\n"
      " scale : float\n"
      "   Real number indicating the scale of the atoms and the cell of the structure. "
          "It is expected, though not enforced, that ``structure.scale * structure.cell`` "
          "is in angstrom and cartesian coordinates. Eventually, such an arbitrary decision "
          "should lead to fewer mistakes.\n"
      " name : string\n"
      "   A string with the name of the structure.\n"  
      " energy : float\n"
      "   Real number indicating the energy of a structure, e.g. for fitting purposes.\n"  
      " weight : float\n"
      "   Real number indicating the weight of a structure in a fitting set, for instance.\n"  
      " freeze : integer\n"
      "   A mask specifying which coordinates are frozen during VFF relaxation.\n"
      " **kwargs\n"
      "   Further keyword arguments are assigned as attributes of the structure instance.\n\n"
      "Note that atoms cannot be set using __init__. "
      "Atoms should be specified after initializing the structure instance, using ``add_atom``.\n",
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

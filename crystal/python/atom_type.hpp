static PyTypeObject AtomStrType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "atom.AtomStr",            /*tp_name*/
    sizeof(AtomStr),             /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)AtomStr_Dealloc, /*tp_dealloc*/
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

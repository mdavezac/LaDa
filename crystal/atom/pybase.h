#if LADA_CRYSTAL_MODULE == 100 || LADA_CRYSTAL_MODULE == 0
  namespace LaDa
  {
    namespace crystal
    {
      namespace // unamed namespace -- declarations unavailable outside of scope.
      {
        //! \brief Describes basic atom type. 
        //! \details This is a python object. 
        struct PyAtomObject
        {
          PyObject_HEAD
          //! Holds list of weak pointers.
          PyObject *weakreflist;
          //! Holds python attribute dictionary.
          PyObject *pydict;
          //! Holds the occupation object.
          PyObject *type;
          //! The atomic position in cartesian coordinate.
          math::rVector3d pos;
        };
#endif 

      
#if LADA_CRYSTAL_MODULE != 1
  //! Returns pointer to atom type.
  PyTypeObject* atom_type()
    LADA_END({ return (PyTypeObject*)api_capsule[LADA_SLOT(crystal)]; })
#else
  api_capsule[LADA_SLOT(crystal)] = (void *)atom_type();
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(LADA_SLOT(crystal))
#include LADA_ASSIGN_SLOT(crystal)
        
#if LADA_CRYSTAL_MODULE != 1
  //! Creates a new atom.
  PyAtomObject* new_atom()
    LADA_END( { return (*(PyAtomObject*(*)())
                        api_capsule[LADA_SLOT(crystal)])(); } )
#else
  api_capsule[LADA_SLOT(crystal)] = (void *)((PyAtomObject*(*)())new_atom);
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(LADA_SLOT(crystal))
#include LADA_ASSIGN_SLOT(crystal)

#if LADA_CRYSTAL_MODULE != 1
  //! Creates a new atom with a given type, also calling initialization.
  PyAtomObject* new_atom(PyTypeObject* _type, PyObject *_args, PyObject *_kwargs)
    LADA_END( { return (*(PyAtomObject*(*)(PyTypeObject*, PyObject*, PyObject*))
                        api_capsule[LADA_SLOT(crystal)])(_type, _args, _kwargs); } )
#else
  api_capsule[LADA_SLOT(crystal)]
         = (void *)((PyAtomObject*(*)(PyTypeObject*, PyObject*, PyObject*))new_atom);
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(LADA_SLOT(crystal))
#include LADA_ASSIGN_SLOT(crystal)

#if LADA_CRYSTAL_MODULE != 1
  //! Creates a deepcopy of atom.
  PyAtomObject *copy_atom(PyAtomObject* _self, PyObject *_memo=NULL)
    LADA_END( { return (*(PyAtomObject*(*)(PyAtomObject*, PyObject*))
                        api_capsule[LADA_SLOT(crystal)])(_self, _memo); })
#else
  api_capsule[3] = (void *)copy_atom;
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(LADA_SLOT(crystal))
#include LADA_ASSIGN_SLOT(crystal)

#if LADA_CRYSTAL_MODULE != 1
  //! Checks type of an object.
  inline bool check_atom(PyObject *_self)
    { return PyObject_TypeCheck(_self, atom_type()); }
  //! Checks type of an object.
  inline bool checkexact_atom(PyObject *_self)
    { return _self->ob_type ==  atom_type(); }
  //! Checks type of an object.
  inline bool check_atom(PyAtomObject *_self)
    { return PyObject_TypeCheck(_self, atom_type()); }
  //! Checks type of an object.
  inline bool checkexact_atom(PyAtomObject *_self)
    { return _self->ob_type ==  atom_type(); }
#endif

#if LADA_CRYSTAL_MODULE != 1
      } // static namespace
    }
  }
#endif

#if LADA_PYTHON_MODULE != 1
  //! Checks whether this a quantity object.
  LADA_INLINE bool check_quantity(PyObject *_in)
    LADA_END( return ((bool(*)(PyObject*))api_capsule[LADA_SLOT(python)])(_in); )
#else
  api_capsule[LADA_SLOT(python)] = (void *)check_quantity;
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(LADA_SLOT(python))
#include LADA_ASSIGN_SLOT(python)


#if LADA_PYTHON_MODULE != 1
  //! Returns a quantity object from a number and a unit in string.
  LADA_INLINE PyObject* fromC_quantity(types::t_real const& _double, std::string const &_units)
    LADA_END( return ( (PyObject*(*)(types::t_real const&, std::string const&) )
                       api_capsule[LADA_SLOT(python)] )
                     (_double, _units); )
#else
  api_capsule[LADA_SLOT(python)] = (void *)
     ( ( PyObject*(*)(types::t_real const &, std::string const&) )  
       fromC_quantity );
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(LADA_SLOT(python))
#include LADA_ASSIGN_SLOT(python)

#if LADA_PYTHON_MODULE != 1
  //! Returns a quantity object from a number and a quantity.
  LADA_INLINE PyObject* fromC_quantity(types::t_real const &_double, PyObject *_unittemplate)
    LADA_END( return ( (PyObject*(*)(types::t_real const &, PyObject*) )  
                       api_capsule[LADA_SLOT(python)] )
                     (_double, _unittemplate); )
#else
  api_capsule[LADA_SLOT(python)] = (void *)
     ( ( PyObject*(*)(types::t_real const &, PyObject*) )  
       fromC_quantity );
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(LADA_SLOT(python))
#include LADA_ASSIGN_SLOT(python)

#if LADA_PYTHON_MODULE != 1
  //! \brief Returns a quantity object from a number and a quantity object.
  //! \details if _number has itself a unit, then it should be convertible to
  //! _units.  However, in that case, _number is returned as is
  //! (Py_INCREF'ed).
  LADA_INLINE PyObject* fromPy_quantity(PyObject *_number, PyObject *_units)
    LADA_END( return ( (PyObject*(*)(PyObject*, PyObject*) )  
                       api_capsule[LADA_SLOT(python)] )
                     (_number, _units); )
#else
  api_capsule[LADA_SLOT(python)] = (void *) fromPy_quantity;
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(LADA_SLOT(python))
#include LADA_ASSIGN_SLOT(python)

#if LADA_PYTHON_MODULE != 1
  //! \brief Returns as a number in specified units.
  //! \details If the number is 0, then one should check whether exception
  //! was thrown.
  LADA_INLINE types::t_real get_quantity(PyObject *_number, std::string const &_units)
    LADA_END( return ( (types::t_real(*)(PyObject*, std::string const&) )  
                       api_capsule[LADA_SLOT(python)] )
                     (_number, _units); )
#else
  api_capsule[LADA_SLOT(python)] = (void *)
     ( ( types::t_real(*)(PyObject*, std::string const&) )  
       get_quantity );
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(LADA_SLOT(python))
#include LADA_ASSIGN_SLOT(python)

#if LADA_PYTHON_MODULE != 1
  //! \brief Returns as a number in specified units.
  //! \details If the number is 0, then one should check whether exception
  //! was thrown.
  LADA_INLINE types::t_real get_quantity(PyObject *_number, PyObject *_units)
    LADA_END( return ( (types::t_real(*)(PyObject*, PyObject*) )  
                       api_capsule[LADA_SLOT(python)] )
                     (_number, _units); )
#else
  api_capsule[LADA_SLOT(python)] = (void *)
     ( ( types::t_real(*)(PyObject*, PyObject*) )  
       get_quantity );
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(LADA_SLOT(python))
#include LADA_ASSIGN_SLOT(python)

#if LADA_PYTHON_MODULE != 1
  //! \brief Returns as float without conversion. 
  LADA_INLINE types::t_real get_quantity(PyObject *_in)
    LADA_END( return ( (types::t_real(*)(PyObject*) )  
                       api_capsule[LADA_SLOT(python)] )
                     (_in); )
#else
  api_capsule[LADA_SLOT(python)] = (void *)
     ( ( types::t_real(*)(PyObject*) ) get_quantity );
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(LADA_SLOT(python))
#include LADA_ASSIGN_SLOT(python)

#ifndef LADA_CRYSTAL_PYTHON_OBJECT_H
#define LADA_CRYSTAL_PYTHON_OBJECT_H

#if LADA_CRYSTAL_MODULE != 1
  class Object;
  //! Object reset function.
  //! Declared as friend to object so that it can be linked at runtime.
  LADA_INLINE void object_reset(PyObject*& _object, PyObject *_in = NULL)
    LADA_END( { return ( *(void(*)(PyObject*&, PyObject*))
                         api_capsule[BOOST_PP_SLOT(1)])(_object, _in); } ) 
#else
  api_capsule[BOOST_PP_SLOT(1)] = (void *)object_reset;
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(BOOST_PP_SLOT(1))
#include BOOST_PP_ASSIGN_SLOT(1)
  
#if LADA_CRYSTAL_MODULE != 1
  LADA_INLINE bool object_equality_op(Object const& _self, Object const &_b)
    LADA_END( { return ( *(bool(*)(Object const&, Object const&))
                         api_capsule[BOOST_PP_SLOT(1)])(_self, _b); } )
#else
  api_capsule[BOOST_PP_SLOT(1)] = (void *)object_equality_op;
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(BOOST_PP_SLOT(1))
#include BOOST_PP_ASSIGN_SLOT(1)

#if LADA_CRYSTAL_MODULE != 1
  //! \brief Dumps representation of an object.
  //! \details Will throw c++ exceptions if python calls fail. Does not clear
  //!          python exceptions.
  LADA_INLINE  std::ostream& operator<< (std::ostream &stream, Object const &_ob)
    LADA_END( { return ( *(std::ostream&(*)(std::ostream&, Object const&))
                         api_capsule[BOOST_PP_SLOT(1)])(_stream, _ob); } )
#else
  api_capsule[BOOST_PP_SLOT(1)]
      = (void *) ( ( (std::ostream&(*)(std::ostream&, Object const&))
                     object_equality_op ) );
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(BOOST_PP_SLOT(1))
#include BOOST_PP_ASSIGN_SLOT(1)

#if LADA_CRYSTAL_MODULE != 1
  //! \brief Thin wrapper around a python refence.
  //! \details In general, steals a reference which it decref's on destruction, unless it
  //!          is null. The destructor is not virtual, hence it is safer to
  //!          keep members of all derived class non-virtual as well.
  //!          When creating this object, the second argument can be false,
  //!          in which case the reference is not owned and never increfed or
  //!          decrefed (eg, borrowed). 
  class Object 
  {
    friend void object_reset(PyObject*& _object, PyObject *_in);
    friend bool object_equality_op(Object const& _self, Object const &_b);
    public:
      //! Steals a python reference.
      Object(PyObject* _in = NULL) : object_(_in) {}
      //! Copies a python reference.
      Object(Object const &_in) { Py_XINCREF(_in.object_); object_ = _in.object_; }
      //! Decrefs a python reference on destruction.
      ~Object() { PyObject* dummy = object_; object_ = NULL; Py_XDECREF(dummy);}
      //! Assignment operator. Gotta let go of current object.
      void operator=(Object const &_in) { Object::reset(_in.borrowed()); }
      //! Assignment operator. Gotta let go of current object.
      void operator=(PyObject *_in) { Object::reset(_in); }
      //! Casts to bool to check validity of reference.
      operator bool() const { return object_ != NULL; }
      //! True if reference is valid.
      bool is_valid () const { return object_ != NULL; }
      //! \brief Resets the wrapped refence.
      //! \details Decrefs the current object if needed. Incref's the input object.
      void reset(PyObject *_in = NULL) { object_reset(_self.object_, _in); }
      //! \brief Resets the wrapped refence.
      //! \details Decrefs the current object if needed.
      void reset(Object const &_in) { reset(_in.object_); }
      
      //! \brief Releases an the reference.
      //! \details After this call, the reference is not owned by this object
      //!          anymore. The reference should be stolen by the caller.
      PyObject* release() { PyObject *result(object_); object_ = NULL; return result; }
      //! Returns a new reference to object.
      PyObject* new_ref() const
      { 
        if(object_ == NULL) return NULL; 
        Py_INCREF(object_);
        return object_; 
      }
      //! Returns a borrowed reference to object.
      PyObject* borrowed() const { return object_; }

      //! \brief Acquires a new reference.
      //! \details First incref's the reference (unless null).
      static Object acquire(PyObject *_in) { Py_XINCREF(_in); return Object(_in); }

      bool hasattr(std::string const &_name) const
        { return PyObject_HasAttrString(object_, _name.c_str()); }

      //! \brief Returns a new reference to a python attribute. 
      //! \details A python exception is set if attribute does not exist, and the
      //!          function returns null.
      inline Object pyattr(std::string const &_name) const
        { return Object(PyObject_GetAttrString((PyObject*)object_, _name.c_str())); }
      //! \brief Returns a new reference to a python attribute. 
      //! \details A python exception is set if attribute does not exist, and the
      //!          function returns null.
      inline Object pyattr(PyObject* _name) const
        { return Object(PyObject_GetAttr((PyObject*)object_, _name)); }
      //! \brief Sets/Deletes attribute.
      inline bool pyattr(std::string const& _name, PyObject *_in)
        { return PyObject_SetAttrString(object_, _name.c_str(), _in) == 0; }
      //! \brief Sets/Deletes attribute.
      inline bool pyattr(std::string const& _name, python::Object const &_in)
        { return PyObject_SetAttrString(object_, _name.c_str(), _in.borrowed()) == 0; }
      //! \brief Sets/Deletes attribute.
      inline bool pyattr(PyObject* _name, PyObject *_in)
        { return PyObject_SetAttr(object_, _name, _in) == 0; }

      //! \brief Compares two objects for equality.
      //! \details Looks for __eq__ in a. If not found, throws c++
      //! exception.
      bool operator==(PyObject *_b) const { return operator==(acquire(_b)); }
      //! \brief Compares two objects for equality.
      //! \details Looks for __eq__ in a. If not found, throws c++
      //! exception.
      bool operator==(Object const &_b) const
        { return object_equality_op(*this, _b); }
      //! \brief Compares two objects for equality.
      //! \details Looks for __eq__ in a. If not found, throws c++
      //! exception.
      bool operator!=(Object const &_b) const { return not operator==(_b); }
      //! \brief Compares two objects for equality.
      //! \details Looks for __eq__ in a. If not found, throws c++
      //! exception.
      bool operator!=(PyObject *_b) const { return operator!=(acquire(_b)); }
    protected:
      //! Python reference.
      PyObject* object_;
  };

  //! \brief Acquires a reference to an object.
  //! \details Input is XINCREF'ed before the return wrappers is created. 
  inline Object acquire(PyObject *_in) { Py_XINCREF(_in); return Object(_in); }
  //! \brief Steals a reference to an object.
  //! \details Input is XINCREF'ed before the return wrappers is created. 
  inline Object steal(PyObject *_in) { return Object(_in); }

  //! \brief Dumps representation of an object.
  //! \details Will throw c++ exceptions if python calls fail. Does not clear
  //!          python exceptions.
  inline std::ostream& operator<< (std::ostream &stream, PyObject* _ob)
    { return stream << Object::acquire(_ob); }
#endif

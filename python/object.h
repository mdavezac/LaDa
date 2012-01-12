#ifndef LADA_PYTHON_OBJECT_H
#define LADA_PYTHON_OBJECT_H

#include <LaDaConfig.h>
#include <Python.h>

namespace LaDa 
{
  namespace python
  {
    //! \brief Thin wrapper around a python refence.
    //! \details In general, steals a reference which it decref's on destruction, unless it
    //!          is null. The destructor is not virtual, hence it is safer to
    //!          keep members of all derived class non-virtual as well.
    //!          When creating this object, the second argument can be false,
    //!          in which case the reference is not owned and never increfed or
    //!          decrefed (eg, borrowed). 
    class Object 
    {
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
        void reset(PyObject *_in = NULL)
        {
          PyObject * const dummy(object_);
          object_ = _in;
          Py_XINCREF(object_);
          Py_XDECREF(dummy);
        }
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

        //! \brief Returns a new reference to a python attribute. 
        //! \details A python exception is set if attribute does not exist, and the
        //!          function returns null.
        inline Object pyattr(std::string const &_name) 
          { return Object(PyObject_GetAttrString((PyObject*)object_, _name.c_str())); }
        //! \brief Returns a new reference to a python attribute. 
        //! \details A python exception is set if attribute does not exist, and the
        //!          function returns null.
        inline Object pyattr(PyObject* _name)
          { return Object(PyObject_GetAttr((PyObject*)object_, _name)); }
        //! \brief   Sets an attribute.
        //! \details If cannot convert object using boost, then returns false
        //!          and sets a python exception.
        template<class T> 
          inline bool pyattr(std::string const &_name, T const &_in) 
          {
            try
            {
              boost::python::object object(_in);
              return PyObject_SetAttrString((PyObject*)object_, _name.c_str(), _in.ptr()) == 0;
            }
            catch(std::exception &e) 
            {
              LADA_PYERROR(internal, ("Could not set atomic attribute " + _name
                                      + ": " + e.what()).c_str() );
              return false;
            }
          }
        //! \brief   Sets an attribute.
        //! \details If cannot convert object using boost, then returns false
        //!          and sets a python exception.
        template<class T> 
          inline bool pyattr(PyObject* _name, T const &_in) 
          {
            try
            {
              boost::python::object object(_in);
              return PyObject_SetAttr((PyObject*)object_, _name, _in.ptr()) == 0;
            }
            catch(std::exception &e) 
            {
              LADA_PYERROR(internal, (std::string("Could not set atomic attribute: ") + e.what()).c_str() );
              return false;
            }
          }
        //! \brief Sets/Deletes attribute.
        inline bool pyattr(std::string const& _name, PyObject *_in)
          { return PyObject_SetAttrString((PyObject*)object_, _name.c_str(), _in) == 0; }
        //! \brief Sets/Deletes attribute.
        inline bool pyattr(PyObject* _name, PyObject *_in)
          { return PyObject_SetAttr((PyObject*)object_, _name, _in) == 0; }
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
    inline std::ostream& operator<< (std::ostream &stream, Object const &_ob)
    {
      PyObject* const repr = PyObject_Repr(_ob.borrowed());
      if(not repr) BOOST_THROW_EXCEPTION(error::internal());
      char const * const result = PyString_AS_STRING(repr);
      if(not result) BOOST_THROW_EXCEPTION(error::internal()); 
      return stream << result;
    }
  }
}
#endif

#ifndef LADA_PYTHON_OBJECT_H
#define LADA_PYTHON_OBJECT_H

#include <LaDaConfig.h>
#include <Python.h>
#include "exceptions.h"

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
        {
          if(not hasattr("__eq__"))
            BOOST_THROW_EXCEPTION( error::TypeError() << 
                                   error::string("No __eq__ member function found "
                                                 "when comparing objects.") );
          Object methodname = PyString_FromString("__eq__");
          if(not methodname)
            BOOST_THROW_EXCEPTION(error::internal() << 
                                  error::string("Could not create string."));
          Object result = PyObject_CallMethodObjArgs(borrowed(), 
                                                     methodname.borrowed(), _b.borrowed(), NULL);
          if(not result) 
            BOOST_THROW_EXCEPTION(error::internal() <<
                                  error::string("Python exception thrown when comparing objects."));
          // Try reflected operation.
          if(result.borrowed() == Py_NotImplemented)
          {
            if(not _b.hasattr("__eq__"))
            {
              if(_b.borrowed()->ob_type != borrowed()->ob_type) return false;
              BOOST_THROW_EXCEPTION( error::TypeError() << 
                                     error::string( "No implementation of "
                                                    "equality between these "
                                                    "two object has been found."));
            }
            result.reset( PyObject_CallMethodObjArgs(_b.borrowed(), 
                                                     methodname.borrowed(), borrowed(), NULL) );
            if(not result)
              BOOST_THROW_EXCEPTION(error::TypeError() << 
                                    error::string("Python exception thrown when comparing objects."));
            if(result.borrowed() == Py_NotImplemented)
            {
              if(_b.borrowed()->ob_type != borrowed()->ob_type) return false;
              BOOST_THROW_EXCEPTION( error::TypeError() <<
                                     error::string( "No implementation of equality between"
                                                    " these two object has been found.") );
            }
          }
          if(PyBool_Check(result.borrowed())) return result.borrowed() == Py_True;
          if(PyInt_Check(result.borrowed())) return PyInt_AS_LONG(result.borrowed()) != 0;
          BOOST_THROW_EXCEPTION( error::ValueError() <<
                                 error::string("Could not make sense of return "
                                               "of comparison function."));
        };
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
    inline std::ostream& operator<< (std::ostream &stream, Object const &_ob)
    {
      PyObject* const repr = PyObject_Repr(_ob.borrowed());
      if(not repr) BOOST_THROW_EXCEPTION(error::internal());
      char const * const result = PyString_AS_STRING(repr);
      if(not result) BOOST_THROW_EXCEPTION(error::internal()); 
      return stream << result;
    }
    //! \brief Dumps representation of an object.
    //! \details Will throw c++ exceptions if python calls fail. Does not clear
    //!          python exceptions.
    inline std::ostream& operator<< (std::ostream &stream, PyObject* _ob)
      { return stream << Object::acquire(_ob); }
  }
}
#endif

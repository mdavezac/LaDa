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
        //! Acquires a python reference.
        Object(PyObject* _in = NULL) : object_(_in) {}
        //! Copies a python reference.
        Object(Object const &_in) : object_(_in.object_) { Py_XINCREF(object_); }
        //! Decrefs a python reference on destruction.
        ~Object() { Py_XDECREF(object_); }
        //! \brief Returns a borrowed reference. 
        //! \details Does not check for null.
        operator PyObject* () const { return object_; }
        //! Casts to bool to check validity of reference.
        operator bool() const { return object_ != NULL; }
        //! True if reference is valid.
        bool is_valid () const { return object_ != NULL; }
        //! \brief Resets the wrapped refence.
        //! \details Decrefs the current object if needed.
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
        PyObject* new_ref() const { Py_XINCREF(object_); return object_; }
        //! Returns a borrowed reference to object.
        PyObject* borrowed() const { return object_; }

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
  }
}
#endif

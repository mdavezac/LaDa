#ifndef LADA_CRYSTAL_ATOM_H
#define LADA_CRYSTAL_ATOM_H

#include "LaDaConfig.h"

#include <boost/python/object.hpp>

#include <python/exceptions.h>

#include "atom_base.h"

namespace LaDa
{
  namespace crystal
  {
    //! \brief Holds a reference to the underlying python atom
    //! \details This wrapper makes it easy to interact with the python object.
    //!          It should be used throughout any c++ code. Mainly, it makes
    //!          sure Py_INCREF and Py_DECREF are called accordingly.
    class Atom
    {
      public: 
        //! Constructor
        Atom() { atom_ = PyAtom_New(); }
        //! Acquires ownership of an atom.
        Atom(Atom const &_c ) { LADA_ACQUIRE_PYOBJECT(atom_, _c.atom_, AtomData); }
        //! Shallow copy Constructor
        Atom(AtomData *_data ) { LADA_ACQUIRE_PYOBJECT(atom_, _data, AtomData); }
        //! Destructor.

        //! Points to owned data.
        AtomData const* operator->() const { return atom_; }
        //! Points to owned data.
        AtomData* operator->() { return atom_; }
        //! Swaps data of two atoms.
        void swap(Atom &_in) { return std::swap(atom_, _in.atom_); }
        //! Acquires another atom. 
        void reset(PyObject *_in) 
        {
          if(_in != NULL and not PyAtom_Check(_in))
          {
            LADA_PYERROR(TypeError, "Cannot acquire object which is not an Atom or subclass.");
            BOOST_THROW_EXCEPTION(error::pyerror() << error::pyobject(_in));
          }
          PyObject *dummy = (PyObject*)atom_;
          atom_ = (AtomData*)_in;
          Py_XINCREF(atom_);
          Py_XDECREF(dummy);
        }
        //! True if the atom is valid.
        bool is_valid() const { return atom_ != NULL; }
        //! \brief Return a (deep) copy of this atom.
        //! \details Return does not share data with this atom. 
        //!          Use constructor to obtain that behavior.
        //!          The user should check that the atom is valid.
        Atom copy() const { return Atom(PyAtom_Copy(atom_)); } 
        //! \brief Returns a new reference to a python attribute. 
        //! \details An python exception is set if attribute does not exist, and the
        //!          function returns null.
        inline PyObject* pyattr(std::string const &_name) 
          { return PyObject_GetAttrString((PyObject*)atom_, _name.c_str()); }
        //! \brief Returns a new reference to a python attribute. 
        //! \details An exception is set if attribute does not exist, and the
        //!          function returns null.
        inline PyObject* pyattr(PyObject* _name)
          { return PyObject_GetAttr((PyObject*)atom_, _name); }
        //! \brief   Sets an attribute.
        //! \details If cannot convert object using boost, then returns false
        //!          and sets a python exception.
        template<class T> 
          inline bool pyattr(std::string const &_name, T const &_in) 
          {
            try
            {
              boost::python::object object(_in);
              return PyObject_SetAttrString((PyObject*)atom_, _name.c_str(), _in.ptr()) == 0;
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
              return PyObject_SetAttr((PyObject*)atom_, _name, _in.ptr()) == 0;
            }
            catch(std::exception &e) 
            {
              LADA_PYERROR(internal, (std::string("Could not set atomic attribute: ") + e.what()).c_str() );
              return false;
            }
          }
        //! \brief Sets/Deletes attribute.
        inline bool pyattr(std::string const& _name, PyObject *_in)
          { return PyObject_SetAttrString((PyObject*)atom_, _name.c_str(), _in) == 0; }
        //! \brief Sets/Deletes attribute.
        inline bool pyattr(PyObject* _name, PyObject *_in)
          { return PyObject_SetAttr((PyObject*)atom_, _name, _in) == 0; }
        //! returns borrowed reference to object.
        PyObject* borrowed() const { return (PyObject*)atom_; }
        //! returns new reference to object.
        PyObject* new_ref() const { Py_XINCREF(atom_); return (PyObject*)atom_; }
    
      private:
        //! Holds a python reference to the atom.
        AtomData* atom_;
    };

    //! \brief Dumps representation of atom.
    //! \details Will throw c++ exceptions if python calls fail. Does not clear
    //!          python exceptions.
    inline std::ostream& operator<< (std::ostream &stream, Atom const &_atom)
    {
      PyObject* const repr = PyObject_Repr(_atom.borrowed());
      if(not repr) BOOST_THROW_EXCEPTION(error::internal());
      char const * const result = PyString_AS_STRING(repr);
      if(not result) BOOST_THROW_EXCEPTION(error::internal()); 
      return stream << result;
    }

  } // namespace Crystal
} // namespace LaDa
  
#endif

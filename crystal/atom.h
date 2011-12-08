#ifndef LADA_CRYSTAL_ATOM_H
#define LADA_CRYSTAL_ATOM_H

#include "LaDaConfig.h"

#include <boost/shared_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/base_object.hpp>

#ifdef LADA_WITH_LNS
#  include "load_n_save/xpr/utilities.h"
#  include "load_n_save/xpr/merge.h"
#endif

#include "atom_base.h"
#ifdef LADA_DO_PYTHON
  #include <Python.h>

  namespace LaDa { namespace crystal { 
    template<class T_TYPE> class Atom;
  } }
  template<class T> PyObject* PyAtom_FromAtom(LaDa::crystal::Atom<T> const &_atom);
#endif
  

namespace LaDa
{
  namespace crystal
  {
    //! \brief Describes an atom where the type is a vector.
    //! \details An atom consists of a position, a type, and frozen status
    //!          variable. The position should always be in cartesian units. The
    //!          type can be anything, from a string with the symbol of the atom,
    //!          to an double wich codes for the atomic type somehow, to a vector
    //!          of strings which describe the possible occupations of the atomic
    //!          position. To this end, the type is a template type \a T_TYPE.
    //!          The frozen status variable indicate which, if any, coordinate
    //!          should not be touched deuring optimization. There are a number
    //!          of possibilities:
    //!            - frozen::None indicate that no coordinate is
    //!                                     frozen.
    //!            - frozen::X indicate that the cartesian x coordinate
    //!                                  is frozen.
    //!            - frozen::Y indicate that the cartesian y coordinate
    //!                                  is frozen.
    //!            - frozen::Z indicate that the cartesian z coordinate
    //!                                  is frozen.
    //!            - frozen::T indicate that the occupation is frozen.
    //!            - Any combination of the above.
    //!            .
    //! \warning The default equality comparison operator compares positions only (not
    //!          occupation).
    template<class T_TYPE>
    class Atom
    {
      friend class boost::serialization::access;
#     ifdef LADA_WITH_LNS
        //! To load and save to xml-like input.
        friend class load_n_save::access; 
#     endif
#     ifdef LADA_DO_PYTHON
        friend PyObject* PyAtom_FromAtom<>(Atom const &_atom);
#     endif
      public: 
        //! Type of the contained data.
        typedef AtomData<T_TYPE> element_type;
        //! Type of the species.
        typedef T_TYPE t_Type;

        //! Constructor
        template<class T_DERIVED>
          Atom   ( Eigen::DenseBase<T_DERIVED> const &_pos, t_Type const &_in,
                   types::t_int _site = -1, types::t_unsigned _freeze = AtomFreezeMixin::frozen::NONE )
               : atom_(new AtomData<T_TYPE>(_pos, _in, _site, _freeze)) {}
        //! Constructor
        Atom() : atom_(new AtomData<T_TYPE>)  {}
        //! Copy Constructor
        Atom(Atom const &_c ) : atom_(_c.atom_) {}
        //! Takes hold of atomic data.
        Atom(boost::shared_ptr<element_type> const& _c) : atom_(_c) {}

        //! Points to owned data.
        AtomData<T_TYPE> const* operator->() const { return atom_.get(); }
        //! Points to owned data.
        AtomData<T_TYPE>* operator->() { return atom_.get(); }
        //! Points to owned data.
        AtomData<T_TYPE>* get() const { return atom_.get(); }
        //! Swaps data of two atoms.
        void swap(Atom &_in) const { return atom_.swap(_in.atom_); }
        //! \brief Return a (deep) copy of this atom.
        //! \details Return does not share data with this atom. 
        //!          Use constructor to obtain that behavior.
        Atom copy() const
          { return Atom<T_TYPE>(atom_->pos, atom_->type, atom_->site, atom_->freeze); }
        //! Returns type.
        T_TYPE const & type() const { return atom_->type; }
        //! Returns type.
        T_TYPE & type() { return atom_->type; }
        //! Returns position.
        math::rVector3d const & pos() const { return atom_->pos; }
        //! Returns position.
        math::rVector3d & pos() { return atom_->pos; }
        //! Returns atomic site index w.r.t to another structure.
        types::t_int const & site() const { return atom_->site; }
        //! Returns atomic site index w.r.t to another structure.
        types::t_int & site() { return atom_->site; }
        //! Returns freeze parameter.
        types::t_unsigned const & freeze() const { return atom_->freeze; }
        //! Returns freeze parameter.
        types::t_unsigned & freeze() { return atom_->freeze; }
#       ifdef LADA_DO_PYTHON
          //! \brief Returns borrowed reference to python dictionary. 
          //! \brief The python dictionary may created at this point, despite
          //!        the const attribute.  Will not throw, but may return NULL
          //!        on error, with a python exception set.
          PyObject* pydict() const
          {
            if(atom_->pydict == NULL) atom_->pydict = PyDict_New();
            return atom_->pydict;
          }
          //! \brief Sets the python attribute dictionary.
          //! \note The reference to _dict is stolen. Your job to increment the
          //!       reference count correctly. The current dictionary, however,
          //!       is decre'f if it exists.
          void pydict(PyObject *_dict) const
          {
            PyObject *dummy = atom_->pydict;
            atom_->pydict = _dict;
            Py_XDECREF(dummy);
          }
          //! Returns a refence to a wrapper around this object.
          PyObject* pyself() const { return PyAtom_FromAtom(*this); }
#       endif
    
      private:
        //! Serializes an atom.
        template<class ARCHIVE> void serialize(ARCHIVE & _ar, const unsigned int _version)
          { _ar & atom_; }
#       ifdef LADA_WITH_LNS
          //! To load and save to xml-like input.
          template<class T_ARCHIVE> bool lns_access(T_ARCHIVE &_ar, const unsigned int _version);
#       endif
        //! \brief This object owns it own data. 
        //! \details Makes python memory management much easier.
        boost::shared_ptr<element_type> atom_;
    };

#   ifdef LADA_WITH_LNS
      //! To load and save to xml-like input.
      template<class T_TYPE> template<class T_ARCHIVE>
        bool Atom<T_TYPE> :: lns_access(T_ARCHIVE &_ar, const unsigned int _version) 
        {
          if(_ar.is_loading() or not atom_) atom_.reset(new AtomData<T_TYPE>());
          return _ar & load_n_save::ext(atom_);
        }
#   endif

    template<class T_TYPE> 
      inline std::ostream& operator<< (std::ostream &stream, Atom<T_TYPE> const &_atom)
        { return stream << *_atom.get(); }

  } // namespace Crystal
} // namespace LaDa
  
#endif

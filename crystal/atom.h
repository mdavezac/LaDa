#ifndef LADA_CRYSTAL_ATOM_H
#define LADA_CRYSTAL_ATOM_H

#include "LaDaConfig.h"

#include <boost/shared_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/base_object.hpp>
#ifdef LADA_DO_PYTHON
#  include <boost/python/object.hpp>
#  include <boost/python/borrowed.hpp>
#endif

#ifdef LADA_WITH_LNS
#  include "load_n_save/xpr/utilities.h"
#  include "load_n_save/xpr/merge.h"
#endif

#include "atom_base.h"

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
      public: 
        //! Type of the species.
        typedef T_TYPE t_Type;

        // Macro to initialize self in constructor.
#       ifdef LADA_SELF
#         error LADA_SELF already defined
#       endif
#       ifdef LADA_DO_PYTHON
#         define LADA_SELF(object) , self_(object)
#       else
#         define LADA_SELF(object)
#       endif
        //! Constructor
        template<class T_DERIVED>
          Atom   ( Eigen::DenseBase<T_DERIVED> const &_pos, t_Type const &_in,
                   types::t_int _site = -1, types::t_unsigned _freeze = AtomFreezeMixin::frozen::NONE )
               : atom_(new AtomData<T_TYPE>(_pos, _in, _site, _freeze)) LADA_SELF(Py_None) {}
        //! Constructor
        Atom() : atom_(new AtomData<T_TYPE>) LADA_SELF(Py_None) {}
        //! Copy Constructor
        Atom(const Atom &_c ) : atom_(_c.atom_) LADA_SELF(_c.self_)
        {
#         ifdef LADA_DO_PYTHON
            if(self_ != Py_None) Py_INCREF(self_)
#         endif
        };
#       undef LADA_SELF
#       ifdef LADA_DO_PYTHON
          //! Constructor
          Atom   (PyObject *_self)
               : atom_(new AtomData<T_TYPE>), self_(_self) {}
            { if(self_ != Py_None) Py_INCREF(self_); };
          //! Returns self object.
          boost::python::object self() const
          { 
            namespace bp = boost::python;
            return bp::object(bp::borrowed<>(self_)); 
          }
          //! Sets python back reference if it is currently null.
          void set_self(PyObject *_self)
          {
            if(self_ != Py_None) return;
            self_ = _self; Py_INCREF(_self); 
          }
#       endif
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
    
      private:
        //! Serializes an atom.
        template<class ARCHIVE> void serialize(ARCHIVE & _ar, const unsigned int _version)
          { _ar & atom_; }
#       ifdef LADA_WITH_LNS
          //! To load and save to xml-like input.
          template<class T_ARCHIVE> bool lns_access(T_ARCHIVE &_ar, const unsigned int _version);
#       endif
#       ifdef LADA_DO_PYTHON
          //! Pointer to python object referencing this object.
          PyObject *self_;
#       endif
          //! \brief This object owns it own data. 
          //! \details Makes python memory management much easier.
          boost::shared_ptr< AtomData<T_TYPE> > atom_;
    };

#   ifdef LADA_WITH_LNS
      //! To load and save to xml-like input.
      template<class T_TYPE> template<class T_ARCHIVE>
        bool Atom<T_TYPE> :: lns_access(T_ARCHIVE &_ar, const unsigned int _version) 
        {
          if(not atom_) atom_.reset(new AtomData<T_TYPE>());
          return _ar & load_n_save::ext(atom_);
        }
#   endif

    template<class T_TYPE> 
      inline std::ostream& operator<< (std::ostream &stream, Atom<T_TYPE> const &_atom)
        { return stream << *_atom.get(); }

  } // namespace Crystal
} // namespace LaDa
  
#endif

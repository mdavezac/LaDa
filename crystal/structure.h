#ifndef LADA_CRYSTAL_STRUCTURE_H
#define LADA_CRYSTAL_STRUCTURE_H

#include "LaDaConfig.h"

#include <boost/shared_ptr.hpp>

#ifdef LADA_WITH_LNS
# include <load_n_save/xpr/utilities.h>
#endif

#include "structure_data.h"


namespace LaDa 
{
  namespace crystal
  {
    // Dummy declaration.
    template<class T_TYPE> class TemplateStructure;

    //! Dumps structure to string.
    template< class T_TYPE >
      std::ostream& operator<<(std::ostream &_stream, TemplateStructure<T_TYPE> const &_str);

    //! Wraps around a shared pointer containing structure data.
    template<class T_TYPE> class TemplateStructure :
       public details::call_add_atom2<TemplateStructure, T_TYPE>
    {
        friend class details::call_add_atom2<TemplateStructure, T_TYPE>;
        friend class boost::serialization::access;
#       ifdef LADA_WITH_LNS
          friend class load_n_save::access;
#       endif
        template<class T> friend
          std::ostream& operator<<(std::ostream &_stream, TemplateStructure<T> const &_str);
        //! Copy Constructor
        TemplateStructure(StructureData<T_TYPE> &_c) : impl_(new StructureData<T_TYPE>(_c)) {}
      public:
        //! \typedef Type of the species
        typedef typename traits::StructureData<T_TYPE>::t_Type t_Type;
        //! \typedef The type of the collection of atoms. 
        typedef typename traits::StructureData<T_TYPE>::t_Atom t_Atom;
        //! \typedef The type of the collection of atoms. 
        typedef typename traits::StructureData<T_TYPE>::t_Atoms t_Atoms;
        //! Type of the iterator.
        typedef typename t_Atoms::iterator iterator;
        //! Type of the constant iterator.
        typedef typename t_Atoms::const_iterator const_iterator;
        //! Type of the reverse iterator.
        typedef typename t_Atoms::reverse_iterator reverse_iterator;
        //! Type of the reverse constant iterator.
        typedef typename t_Atoms::const_reverse_iterator const_reverse_iterator;
        //! Type of the atoms.
        typedef typename t_Atoms::value_type value_type;
        //! Type of the reference to the atoms.
        typedef typename t_Atoms::reference reference;
        //! Type of the constant reference to the atoms.
        typedef typename t_Atoms::const_reference const_reference;
        //! Type of the size.
        typedef typename t_Atoms::size_type size_type;
        //! Type of the differnce.
        typedef typename t_Atoms::difference_type difference_type;
        //! Type of the pointer to the atoms.
        typedef typename t_Atoms::pointer pointer;
        //! Type of the allocator used for the atoms.
        typedef typename t_Atoms::allocator_type allocator_type;

        //! Constructor
        TemplateStructure() : impl_(new StructureData<T_TYPE>()) {};
        //! Cloning Constructor
        TemplateStructure(const TemplateStructure &_c) : impl_(_c.impl_) {}
        //! Destructor.
        ~TemplateStructure () {};

        //! Returns const reference to cell.
        math::rMatrix3d const & cell() const { return impl_->cell; }
        //! Returns reference to cell.
        math::rMatrix3d & cell() { return impl_->cell; }
        //! Returns const reference to name.
        std::string const & name() const { return impl_->name; }
        //! Returns reference to name.
        std::string & name() { return impl_->name; }
        //! Returns const reference to energy.
        types::t_real const & energy() const { return impl_->energy; }
        //! Returns reference to energy.
        types::t_real & energy() { return impl_->energy; }
        //! Returns const reference to weight.
        types::t_real const & weight() const { return impl_->weight; }
        //! Returns reference to weight.
        types::t_real & weight() { return impl_->weight; }
        //! Returns const reference to scale.
        types::t_real const & scale() const { return impl_->scale; }
        //! Returns reference to scale.
        types::t_real & scale() { return impl_->scale; }
        //! Returns const reference to freeze (frozen degrees of freedom).
        types::t_unsigned const & freeze() const { return impl_->freeze; }
        //! Returns reference to freeze.
        types::t_unsigned & freeze() { return impl_->freeze; }

        //! Deep copy of a structure.
        TemplateStructure copy() const { return TemplateStructure(*impl_); }
        //! Swaps content of two structures.
        void swap(TemplateStructure &_other) { impl_.swap(_other.impl_); }

        //! Iterator to the atoms.
        iterator begin() { return impl_->atoms.begin(); }
        //! Iterator to the atoms.
        iterator end() { return impl_->atoms.end(); }
        //! Iterator to the atoms.
        reverse_iterator rbegin() { return impl_->atoms.rbegin(); }
        //! Iterator to the atoms.
        reverse_iterator rend() { return impl_->atoms.rend(); }
        //! Iterator to the atoms.
        const_iterator begin() const { return impl_->atoms.begin(); }
        //! Iterator to the atoms.
        const_iterator end() const { return impl_->atoms.end(); }
        //! Iterator to the atoms.
        const_reverse_iterator rbegin() const { return impl_->atoms.rbegin(); }
        //! Iterator to the atoms.
        const_reverse_iterator rend() const { return impl_->atoms.rend(); }
        //! Number of atoms.
        size_type size() const { return impl_->atoms.size(); }
        //! Maximum number of atoms.
        size_type max_size() const { return impl_->atoms.max_size(); }
        //! Number of atoms.
        void resize(size_type _n) { impl_->atoms.size(_n); }
        //! Number of atoms.
        void resize(size_type _n, const_reference _obj) { impl_->atoms.size(_n, _obj); }
        //! Attempts to reserve memory for atoms.
        void reserve(size_type _n) { impl_->atoms.reserve(_n); }
        //! Reserved memory.
        size_type capacity() const { return impl_->atoms.capacity(); }
        //! Whether any atoms are in the structure.
        bool empty() const { return impl_->empty(); }
        //! Returns nth atom.
        reference operator[](size_type _n) { return impl_->atoms[_n]; }
        //! Returns nth atom.
        const_reference operator[](size_type _n) const { return impl_->atoms[_n]; }
        //! Returns nth atom.
        reference at(size_type _n) { return impl_->atoms.at(_n); }
        //! Returns nth atom.
        const_reference at(size_type _n) const { return impl_->atoms.at(_n); }
        //! Returns nth atom.
        reference front() { return impl_->atoms.front(); }
        //! Returns nth atom.
        const_reference front() const { return impl_->atoms.front(); }
        //! Returns nth atom.
        reference back() { return impl_->atoms.back(); }
        //! Returns nth atom.
        const_reference back() const { return impl_->atoms.back(); }
        //! Replaces content of the container.
        template <class InputIterator>
          void assign(InputIterator _first, InputIterator _last)
          { impl_->atoms.assign(_first, _last); }
        //! Replaces content of the container.
        void assign(size_type _n, const_reference _u) { impl_->atoms.assign(_n, _u); }
        //! Adds atom at end of container.
        void push_back(const_reference _u) { impl_->atoms.push_back(_u); }
        //! Deletes last element. 
        void pop_back() { impl_->atoms.pop_back(); }
        //! Inserts atoms in container at given position.
        template <class InputIterator>
          void insert(iterator _pos, InputIterator _first, InputIterator _last)
          { impl_->atoms.insert(_pos, _first, _last); }
        //! Inserts one atom at given position.
        iterator insert(iterator _pos, const_reference x) { impl_->atoms.insert(_pos, x); }
        //! Inserts n atom at given position.
        void insert(iterator _pos, size_type _n, const_reference x)
          { impl_->atoms.insert(_pos, _n, x); }
        //! Erases atom at given positions.
        iterator erase(iterator _pos) { return impl_->atoms.erase(_pos); }
        //! Erases range of atoms at given positions.
        iterator erase(iterator _first, iterator _last) { return impl_->atoms.erase(_first, _last); }
        //! Clears all atoms from structure.
        void clear() { impl_->atoms.clear(); }
        //! Returns allocator object used to construct atom container.
        allocator_type get_allocator() const { return impl_->atoms.get_allocator(); }
        //! Initializer for cell.
        math::details::SetCell< boost::mpl::int_<1> >
          set_cell(types::t_real _x, types::t_real _y, types::t_real _z)
            { return impl_->set_cell(_x, _y, _z); }
        //! Initializer for cell.
        math::details::SetCell< boost::mpl::int_<1> >
          set_cell(math::rVector3d _pos)
            { return impl_->set_cell(_pos); }

        //! Access to cell parameters
        types::t_real operator()(size_t i, size_t j) const { return cell()(i,j); }
        //! Access to cell parameters
        types::t_real& operator()(size_t i, size_t j) { return cell()(i,j); }

        //! \brief True if both structures refer to the same object in memory.
        //! \details Does not compare values, just memory objects.
        bool is_same(TemplateStructure const &_in) { return impl_ == _in.impl; }

        //! Returns  structure volume.
        types::t_real volume() const { return std::abs(impl_->cell.determinant()); }

        //! Transforms a structure according to an affine transformation.
        TemplateStructure transform(math::Affine3d const &_affine) const;
      private:
        //! Serializes a structure.
        template<class ARCHIVE> void serialize(ARCHIVE & _ar, const unsigned int _version)
          { _ar & impl_; }
#       ifdef LADA_WITH_LNS
          //! To load and save to xml-like input.
          template<class T_ARCHIVE>
            bool lns_access(T_ARCHIVE &_ar, load_n_save::version_type const _version);
#       endif
        //! Holds data.
        boost::shared_ptr< StructureData<T_TYPE> > impl_;
    };

    template< class T_TYPE >
      std::ostream& operator<<(std::ostream &_stream, TemplateStructure<T_TYPE> const &_str)
        { return _stream << *_str.impl_; }

#   ifdef LADA_WITH_LNS
      template<class T_TYPE> template<class T_ARCHIVE>
        bool TemplateStructure<T_TYPE> :: lns_access(T_ARCHIVE &_ar, load_n_save::version_type const _version) 
        {
          if(not impl_) impl_.reset(new StructureData<T_TYPE>());
          return _ar & load_n_save::ext(impl_);
        }
#   endif
    //! Transforms a structure according to an affine transformation.
    template<class T_TYPE> 
      TemplateStructure<T_TYPE> TemplateStructure<T_TYPE>::transform(math::Affine3d const &_affine) const
      {
        TemplateStructure<T_TYPE> result = copy();
        result.cell() = _affine.linear() * cell();
        iterator i_first = result.begin();
        iterator const i_end = result.end();
        for(; i_first != i_end;  ++i_first)
          i_first->pos = _affine * i_first->pos;
        return result;
      }

  } // namespace crystal
} // namespace LaDa

#endif

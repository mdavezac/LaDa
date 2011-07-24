#ifndef LADA_CRYSTAL_STRUCTURE_H
#define LADA_CRYSTAL_STRUCTURE_H

#include "LaDaConfig.h"

#include "structure_data.h"


namespace LaDa 
{
  namespace crystal
  {
    // Dummy declaration.
    template<class TYPE> class TemplateStructure;

    //! Dumps structure to string.
    template< class T_TYPE >
      std::ostream& operator<<(std::ostream &_stream, TemplateStructure<T_TYPE> const &_str);

    //! Wraps around a shared pointer containing structure data.
    template<class TYPE> class TemplateStructure 
    {
        friend class boost::serialization::access;
#       ifdef LADA_WITH_LNS
          friend class load_n_save::access;
#       endif
        friend template< class T_TYPE >
          std::ostream& operator<<(std::ostream &_stream, TemplateStructure<T_TYPE> const &_str);
      public:
        //! Namespace for the frozen dof.
        typedef typename StructureData<TYPE> :: frozen frozen;
        //! Type of the iterator.
        typedef typename StructureData<TYPE>::t_Atoms::iterator iterator;
        //! Type of the constant iterator.
        typedef typename StructureData<TYPE>::t_Atoms::const_iterator const_iterator;
        //! Type of the reverse iterator.
        typedef typename StructureData<TYPE>::t_Atoms::reverse_iterator reverse_iterator;
        //! Type of the reverse constant iterator.
        typedef typename StructureData<TYPE>::t_Atoms::reverse_const_iterator reverse_const_iterator;
        //! Type of the atoms.
        typedef typename StructureData<TYPE>::t_Atoms::value_type value_type;
        //! Type of the reference to the atoms.
        typedef typename StructureData<TYPE>::t_Atoms::reference reference;
        //! Type of the constant reference to the atoms.
        typedef typename StructureData<TYPE>::t_Atoms::const_reference const_reference;
        //! Type of the size.
        typedef typename StructureData<TYPE>::t_Atoms::size_type size_type;
        //! Type of the differnce.
        typedef typename StructureData<TYPE>::t_Atoms::difference_type difference_type;
        //! Type of the pointer to the atoms.
        typedef typename StructureData<TYPE>::t_Atoms::pointer_type pointer_type;
        //! Type of the constant pointer to the atoms.
        typedef typename StructureData<TYPE>::t_Atoms::const_pointer_type const_pointer_type;
        //! Type of the allocator used for the atoms.
        typedef typename StructureData<TYPE>::t_Atoms::allocator_type allocator_type;

        //! Constructor
        TemplateStructure() : impl_(new StructureData<TYPE>()) {};
        //! Copy Constructor
        TemplateStructure   ( const TemplateStructure &_str )
                          : impl_(new StructureData<TYPE>(_c.impl_)) {}
        //! Destructor.
        ~TemplateStructure () {};

        //! Returns const reference to cell.
        rMatrix3d const & cell() const { return impl_->cell; }
        //! Returns reference to cell.
        rMatrix3d & cell() const { return impl_->cell; }
        //! Returns const reference to name.
        std::string const & name() const { return impl_->name; }
        //! Returns reference to name.
        std::string & name() const { return impl_->name; }
        //! Returns const reference to energy.
        types::t_real const & energy() const { return impl_->energy; }
        //! Returns reference to energy.
        types::t_real & energy() const { return impl_->energy; }
        //! Returns const reference to weight.
        types::t_real const & weight() const { return impl_->weight; }
        //! Returns reference to weight.
        types::t_real & weight() const { return impl_->weight; }
        //! Returns const reference to scale.
        types::t_real const & scale() const { return impl_->scale; }
        //! Returns reference to scale.
        types::t_real & scale() const { return impl_->scale; }
        //! Returns const reference to freeze (frozen degrees of freedom).
        types::t_unsigned const & freeze() const { return impl_->freeze; }
        //! Returns reference to freeze.
        types::t_unsigned & freeze() const { return impl_->freeze; }

        //! Shallow copy, e.g. refers to same data.
        TemplateStructure shallow_copy() const
        {
          TemplateStructure result; result.impl_ = impl_; 
          return result;
        }
        //! Swaps content of two structures.
        void swap(TemplateStructure &_structure) { impl_.swap(_other.impl_); }

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
        reverse_const_iterator rbegin() const { return impl_->atoms.rbegin(); }
        //! Iterator to the atoms.
        reverse_const_iterator rend() const { return impl_->atoms.rend(); }
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
        bool empty() const { return impbool empty(); }
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


      private:
        //! Serializes a structure.
        template<class ARCHIVE> void serialize(ARCHIVE & _ar, const unsigned int _version)
          { ar & impl_; }
#     ifdef LADA_WITH_LNS
        //! To load and save to xml-like input.
        template<class T_ARCHIVE>
          bool lns_access(T_ARCHIVE &_ar, load_n_save::version_type const _version);
#     endif
        //! Holds data.
        boost::shared_ptr< StructureData<TYPE> > impl_;
    };

    template< class T_TYPE >
      std::ostream& operator<<(std::ostream &_stream, TemplateStructure<T_TYPE> const &_str)
        { return _stream << *_str.impl_; }

#   ifdef LADA_WITH_LNS
      template<class TYPE> template<class T_ARCHIVE>
        bool TemplateStructure :: lns_access(T_ARCHIVE &_ar, load_n_save::version_type const _version) 
        {
          if( _ar.is_loading() )
          {
            boost::shared_ptr< StructureData<TYPE> > dummy(impl_);
            impl_.reset(new StructureData<TYPE>());
          }
          return _ar & *impl_;
        }
#   endif
  } // namespace crystal
} // namespace LaDa

#endif

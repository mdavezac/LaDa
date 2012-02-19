#ifndef LADA_CRYSTAL_STRUCTURE_H
#define LADA_CRYSTAL_STRUCTURE_H

#include "LaDaConfig.h"

#define PY_ARRAY_UNIQUE_SYMBOL crystal_ARRAY_API
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>


#include "pybase.h"


namespace LaDa 
{
  namespace crystal
  {
    //! Wraps a python structure. 
    class Structure : public python::Object
    {
      public:
        //! \typedef The type of the collection of atoms. 
        typedef Atom t_Atom;
        //! \typedef The type of the collection of atoms. 
        typedef std::vector<t_Atom> t_Atoms;
        //! \typedef Type of the iterator.
        typedef t_Atoms::iterator iterator;
        //! \typedef Type of the constant iterator.
        typedef t_Atoms::const_iterator const_iterator;
        //! \typedef Type of the reverse iterator.
        typedef t_Atoms::reverse_iterator reverse_iterator;
        //! \typedef Type of the reverse constant iterator.
        typedef t_Atoms::const_reverse_iterator const_reverse_iterator;
        //! \typedef Type of the atoms.
        typedef t_Atoms::value_type value_type;
        //! \typedef Type of the reference to the atoms.
        typedef t_Atoms::reference reference;
        //! \typedef Type of the constant reference to the atoms.
        typedef t_Atoms::const_reference const_reference;
        //! \typedef Type of the size.
        typedef t_Atoms::size_type size_type;
        //! \typedef Type of the difference.
        typedef t_Atoms::difference_type difference_type;
        //! \typedef Type of the pointer to the atoms.
        typedef t_Atoms::pointer pointer;
        //! \typedef Type of the allocator used for the atoms.
        typedef t_Atoms::allocator_type allocator_type;

        //! Constructor
        Structure() : Object() { object_ = (PyObject*) PyStructure_New(); }
        //! Cloning Constructor
        Structure(const Structure &_c) : Object(_c) {}
        //! Shallow copy Constructor
        Structure(StructureData *_data ) : Object((PyObject*)_data) {}
        //! Full Initialization.
        Structure(PyObject *_args, PyObject *_kwargs) : Object()
          { object_ = (PyObject*)PyStructure_NewFromArgs(structure_type(), _args, _kwargs); }

        //! Returns const reference to cell.
        math::rMatrix3d const & cell() const { return ((StructureData*)object_)->cell; }
        //! Returns reference to cell.
        math::rMatrix3d & cell() { return ((StructureData*)object_)->cell; }
        //! Returns const reference to cell.
        math::rMatrix3d::Scalar const & cell(size_t i, size_t j) const
          { return ((StructureData*)object_)->cell(i, j); }
        //! Returns reference to cell.
        math::rMatrix3d::Scalar & cell(size_t i, size_t j)
          { return ((StructureData*)object_)->cell(i, j); }
        //! Returns const reference to scale.
        types::t_real const & scale() const { return ((StructureData*)object_)->scale; }
        //! Returns reference to scale.
        types::t_real & scale() { return ((StructureData*)object_)->scale; }

        //! Deep copy of a structure.
        Structure copy() const { return Structure(PyStructure_Copy((StructureData*)object_)); }
        //! Swaps content of two structures.
        void swap(Structure &_other) { std::swap(object_, _other.object_); }
        //! Points to data.
        StructureData const* operator->() const { return (StructureData*)object_; }
        //! Points to data.
        StructureData* operator->() { return (StructureData*)object_; }

        //! Iterator to the atoms.
        iterator begin() { return ((StructureData*)object_)->atoms.begin(); }
        //! Iterator to the atoms.
        iterator end() { return ((StructureData*)object_)->atoms.end(); }
        //! Iterator to the atoms.
        reverse_iterator rbegin() { return ((StructureData*)object_)->atoms.rbegin(); }
        //! Iterator to the atoms.
        reverse_iterator rend() { return ((StructureData*)object_)->atoms.rend(); }
        //! Iterator to the atoms.
        const_iterator begin() const { return ((StructureData*)object_)->atoms.begin(); }
        //! Iterator to the atoms.
        const_iterator end() const { return ((StructureData*)object_)->atoms.end(); }
        //! Iterator to the atoms.
        const_reverse_iterator rbegin() const { return ((StructureData*)object_)->atoms.rbegin(); }
        //! Iterator to the atoms.
        const_reverse_iterator rend() const { return ((StructureData*)object_)->atoms.rend(); }
        //! Number of atoms.
        size_type size() const { return ((StructureData*)object_)->atoms.size(); }
        //! Maximum number of atoms.
        size_type max_size() const { return ((StructureData*)object_)->atoms.max_size(); }
        //! Number of atoms.
        void resize(size_type _n) { ((StructureData*)object_)->atoms.resize(_n); }
        //! Number of atoms.
        void resize(size_type _n, const_reference _obj) { ((StructureData*)object_)->atoms.resize(_n, _obj); }
        //! Attempts to reserve memory for atoms.
        void reserve(size_type _n) { ((StructureData*)object_)->atoms.reserve(_n); }
        //! Reserved memory.
        size_type capacity() const { return ((StructureData*)object_)->atoms.capacity(); }
        //! Whether any atoms are in the structure.
        bool empty() const { return ((StructureData*)object_)->atoms.empty(); }
        //! Returns nth atom.
        reference operator[](size_type _n) { return ((StructureData*)object_)->atoms[_n]; }
        //! Returns nth atom.
        const_reference operator[](size_type _n) const { return ((StructureData*)object_)->atoms[_n]; }
        //! Returns nth atom.
        reference at(size_type _n) { return ((StructureData*)object_)->atoms.at(_n); }
        //! Returns nth atom.
        const_reference at(size_type _n) const { return ((StructureData*)object_)->atoms.at(_n); }
        //! Returns nth atom.
        reference front() { return ((StructureData*)object_)->atoms.front(); }
        //! Returns nth atom.
        const_reference front() const { return ((StructureData*)object_)->atoms.front(); }
        //! Returns nth atom.
        reference back() { return ((StructureData*)object_)->atoms.back(); }
        //! Returns nth atom.
        const_reference back() const { return ((StructureData*)object_)->atoms.back(); }
        //! Replaces content of the container.
        template <class InputIterator>
          void assign(InputIterator _first, InputIterator _last)
          { ((StructureData*)object_)->atoms.assign(_first, _last); }
        //! Replaces content of the container.
        void assign(size_type _n, const_reference _u) { ((StructureData*)object_)->atoms.assign(_n, _u); }
        //! Adds atom at end of container.
        void push_back(const_reference _u) { ((StructureData*)object_)->atoms.push_back(_u); }
        //! Deletes last element. 
        void pop_back() { ((StructureData*)object_)->atoms.pop_back(); }
        //! Inserts atoms in container at given position.
        template <class InputIterator>
          void insert(iterator _pos, InputIterator _first, InputIterator _last)
          { ((StructureData*)object_)->atoms.insert(_pos, _first, _last); }
        //! Inserts one atom at given position.
        iterator insert(iterator _pos, const_reference x) { ((StructureData*)object_)->atoms.insert(_pos, x); }
        //! Inserts n atom at given position.
        void insert(iterator _pos, size_type _n, const_reference x)
          { ((StructureData*)object_)->atoms.insert(_pos, _n, x); }
        //! Erases atom at given positions.
        iterator erase(iterator _pos) { return ((StructureData*)object_)->atoms.erase(_pos); }
        //! Erases range of atoms at given positions.
        iterator erase(iterator _first, iterator _last)
          { return ((StructureData*)object_)->atoms.erase(_first, _last); }
        //! Clears all atoms from structure.
        void clear() { ((StructureData*)object_)->atoms.clear(); }
        //! Returns allocator object used to construct atom container.
        allocator_type get_allocator() const { return ((StructureData*)object_)->atoms.get_allocator(); }

        //! \brief True if both structures refer to the same object in memory.
        //! \details Does not compare values, just memory objects.
        bool is_same(Structure const &_in) { return object_ == _in.object_; }

        //! Returns  structure volume.
        types::t_real volume() const
          { return std::abs(((StructureData*)object_)->cell.determinant()) 
                   * std::pow(((StructureData*)object_)->scale, 3); }

        //! Transforms a structure according to an affine transformation.
        void transform(Eigen::Matrix<types::t_real, 4, 3> const &_affine);

        //! Returns borrowed reference to dictionary.
        PyObject* dict() const { return ((StructureData*)object_)->pydict; }

        //! Check if instance is an atom.
        static bool check(PyObject *_str) { return PyStructure_Check(_str); }
        //! Check if instance is an atom.
        static bool check_exact(PyObject *_str) { return PyStructure_CheckExact(_str); }
        //! \brief Acquires new reference to an object.
        //! \details incref's reference first, unless null.
        //!          If non-null checks that it is a subtype of Structure.
        //! \throws error::TypeError if not an Structure or subtype, both cpp and python.
        static Structure acquire(PyObject *_str) 
        {
          if(_str == NULL) return Structure((StructureData*)_str);
          if(not Structure::check(_str))
          {
            LADA_PYERROR_FORMAT( TypeError,
                                 "Expected an Structure or subtype, not %.200s",
                                 _str->ob_type->tp_name );
            BOOST_THROW_EXCEPTION(error::TypeError());
          }
          Py_INCREF(_str);
          return Structure((StructureData*)_str);
        }
        //! \brief Acquires new reference to an object.
        //! \details incref's reference first, unless null.
        //!          Does not check object is Structure or subtype, and does not
        //!          throw.
        static Structure acquire_(PyObject *_str) 
          { Py_XINCREF(_str); return Structure((StructureData*)_str); }
    };

  } // namespace crystal
} // namespace LaDa

#endif

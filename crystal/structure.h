#ifndef _ISING_CE_STRUCTURE_H_
#define _ISING_CE_STRUCTURE_H_

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
          bool lns_access(T_ARCHIVE &_ar, load_n_save::version_type const _version) 
          {
            if( _ar.is_loading() )
            {
              boost::shared_ptr< StructureData<TYPE> > dummy(impl_);
              impl_.reset(new StructureData<TYPE>());
            }
            return _ar & *impl_;
          }
#     endif
        //! Holds data.
        boost::shared_ptr< StructureData<TYPE> > impl_;
    };

    struct Structure : public TStructure< types::t_real >
    {
      friend class boost::serialization::access;
      //! The atomic type
      typedef Atom_Type<types::t_real>  t_Atom;
      //! The type of the collection of atoms. 
      typedef std::vector< t_Atom >     t_Atoms;
      //! The reciprocal-space vector type
      typedef CAtom                     t_kAtom;
      //! The type of the collection of reciprocal-space vector.
      typedef std::vector< t_kAtom >    t_kAtoms;

      //! The reciprocal-space vector position in cartesian unit and their intensity.
      std::vector< CAtom > k_vecs;

      public: 

      //! Constructor
      Structure() : TStructure<types::t_real>() {}
      //! Constructor. Loads itself from XML node \a _element.
      Structure   ( const TiXmlElement &_element )
                : TStructure<types::t_real>() { Load( _element ); };
      //! Copy Constructor
      Structure   ( const Structure &_str )
                : TStructure<types::t_real>( _str ), k_vecs(_str.k_vecs) {}
      //! Destructor.
      ~Structure () {};


      //! Prints a structure to a string
      const std::string string() const
        { std::ostringstream sstr; print_out(sstr); return sstr.str(); }

      //! \brief Sets the occupations from the vector \a %types.
      //! \details If \a %types is smaller than Structure::atoms than only the
      //!          first occupations are set. If is is larger, than only the
      //!          first components of \a %types are used.
      //! \pre The atoms in this structure and the occupations of \a %types should
      //!      be listed in the same order.
      void set_atom_types( const std::vector<types::t_real> &types);
      //! \brief Copies the atomic occupations ton \a %types.
      //! \details If \a %types is smaller than Structure::atoms, then only the
      //!          first atoms appearing in Structure::atoms are copied. If \a
      //!          %types is larger than  Structure::atoms, then only the first
      //!          components of \a %types are set.
      void get_atom_types( std::vector<types::t_real> &types) const;
      //! Returns the average occupation.
      types::t_real get_concentration() const;
      //! Loads a structure from XML.
      bool Load( const TiXmlElement &_element );
      //! Saves a structure to XML.
      void print_xml( TiXmlElement &_node ) const;

      //! \brief Copies the reciprocal-space vectors to a container.
      //! \details Since Atom_Type automatically returns reference to a
      //!          math::rVector3d and its type, this routien can copy the full
      //!          reciprocal space vector, the positions only, or the
      //!          intensities only, depending on the type of
      //!          t_container::value_type.
      template<class t_container > void set_kvectors( const t_container &_container );

      //! \brief Computes the position of the reciprocal-space vectors in the first
      //!        Brillouin zone of the lattice unit-cell.
      void find_k_vectors();
      //! \brief Prints in XCrysDen format 
      //! \see <A HREF="http://www.xcrysden.org">www.xcrysden.org</A>.
      std::ostream& print_xcrysden( std::ostream &_stream ) const;
      //! \brief Prints in xyz format 
      //! \details XYZ format is as follows
      //! \code
      //     2
      //     Diamond unit-cell
      //     C  0 0 0
      //     C  0.25 0.25 0.25
      //! \endcode
      //! The first line contains the number of atoms.
      //! The second line contains a comment or name.
      //! The following lines define each atom throught its symbol and three
      //! coordinates. There should be as many lines as specified above.
      //! \param[in] _stream where to output the structure in xyz
      //! \param[in] _name second line of xyz format is a name or comment. Fill
      //!                  in here. This string should not include any endline
      //!                  character.
      std::ostream& print_xyz( std::ostream &_stream,
                               const std::string &_name = "" ) const;
      
#     ifdef LADA_WITH_LNS
        using TStructure<types::t_real> :: lns_access;
#     endif

      protected:
        //! Finds parent node.
        const TiXmlElement* find_node( const TiXmlElement &_element );
        //! Loads attributes and non-essential variables.
        void load_attributes( const TiXmlElement &_element );
        //! Loads the cell.
        bool load_cell( const TiXmlElement &_element );
        //! loads atoms.
        bool load_atoms( const TiXmlElement &_element );
        //! loads k-vectors.
        bool load_kvecs( const TiXmlElement &_element );
        //! loads an epitaxial structure.
        bool load_epitaxial( const TiXmlElement &_node );

      private:
        //! Serializes a structure.
        template<class ARCHIVE> void serialize(ARCHIVE & _ar, const unsigned int _version);
    };

    //! Returns true if \a _a and \a _b are periodic equivalents of the unit-cell \a _cell.
    bool are_equivalent( const math::rVector3d &_a,
                         const math::rVector3d &_b,
                         const math::rMatrix3d &_cell);
    //! Dumps a structure to a stream.
    template< class T_TYPE >
      std::ostream& operator<<( std::ostream& _stream, const Crystal::TStructure<T_TYPE>& _struc )
        { _struc.print_out(_stream); return _stream; }

    //! compares two kvectors according to length and position.
    bool sort_kvec( const math::rVector3d &_vec1, const math::rVector3d &_vec2 );

    //! \brief Creates an epitaxial structure.
    //! \param[out] _structure the structure on output. The atoms are ordered
    //!                        with respect to the growth direction.
    //! \param[in] _direction the growth direction.
    //! \param[in] _extent the size in the direction growth ( x ), the other two
    //!                    directions. These are taken from the lattice
    //!                    unit-cell such that the determinant of the structure
    //!                    cell is strictly positive. The first vector in the
    //!                    unit-cell is the growth direction.
    bool create_epitaxial_structure( Structure& _structure,
                                     math::rVector3d &_direction,
                                     math::iVector3d &_extent );

    //! Returns concentration over set site.
    types::t_real concentration( const Structure& _structure, const size_t i );

    //! Returns the tag for freezing coordinate i, j of structure cell.
    inline types::t_unsigned cell_freeze_tag( const size_t i, const size_t j )
    {
      if( i == 0 and j == 0 ) return Structure :: FREEZE_XX;
      if( (i == 0 and j == 1) or (i == 1 and j == 0) ) return Structure :: FREEZE_XY;
      if( (i == 0 and j == 2) or (i == 2 and j == 0) ) return Structure :: FREEZE_XZ;
      if( i == 1 and j == 1 ) return Structure :: FREEZE_YY;
      if( i == 1 and j == 2 ) return Structure :: FREEZE_YZ;
      if( i == 2 and j == 2 ) return Structure :: FREEZE_ZZ;
      return 0;
    };

    void convert_real_to_string_structure( Structure const& _real,
                                           TStructure<std::string> &_string );
    void convert_string_to_real_structure( TStructure<std::string> const &_string,
                                           Structure &_real );

  } // namespace Crystal

} // namespace LaDa

#include "structure.impl.h"


#endif

#ifndef _VFF_ATOMIC_CENTER_H_
#define _VFF_ATOMIC_CENTER_H_

#include "LaDaConfig.h"

#include <vector>

#include <tinyxml/tinyxml.h>

#include <crystal/atom.h>
#include <crystal/structure.h>
#include <crystal/lattice.h>
#include <opt/types.h>
#include <opt/function_base.h>
#include <opt/debug.h>

#include "data.h"

namespace LaDa
{
  namespace vff
  {

    class Vff;

    //! \brief Represents a single Structure::t_Atom and its first neighbor relationships
    //! \details This class is meant to be used in conjunction with a list of
    //! Crystal::TStructure<std::string>::t_Atom, most likely in a Crystal::TStructure<std::string>. It contains a
    //! pointer, AtomicCenter::i_atom_, which points a single Crystal::TStructure<std::string>::t_Atom. The
    //! first neighbor bonds of this atom are collected as vector of pointers to
    //! Vff::AtomicCenter objects in AtomicCenter::bonds. Since we are concerned
    //! with periodic structures, AtomicCenter::translations and
    //! AtomicCenter::bool record which periodic image of an atom
    //! AtomicCenter::bonds refer to.
    class AtomicCenter
    {
      friend class Vff;
      //! The type of the atom  
      typedef Crystal::TStructure<std::string>::t_Atom  t_Atom;
      //! The container of atomic centers. Defined here once and for all.
      typedef std::vector<AtomicCenter> t_Centers;
      //! Type of pointer/iterator to the atomic center on the other side of the bond
      typedef t_Centers :: iterator t_Bond;
      //! A reference to the of pointer/iterator to the atomic center on the other side of the bond
      typedef t_Centers :: iterator& t_BondRefd;
      //! \brief A constant reference to the of pointer/iterator to the atomic
      //!        center on the other side of the bond
      typedef const t_Centers :: iterator& const_t_BondRefd;
      //! A transformation from a t_Center :: iterator to a t_Bond
      static t_BondRefd __make__iterator__( t_Centers::iterator &_i ) { return _i; }

      public:
        class const_iterator;
        
      protected:
        t_Atom *i_atom_; //!< The atom this object is addressing
        //! \brief Other Vff::AtomicCenter objects with which this one is in a bond-relationship
        //! \details Via bonds, a collection of AtomicCenter can be made into a tree,
        //! which can be travelled linearly, or through the first neighbor bonds,
        //! using AtomicCenter::const_iterator.
        //! \sa AtomicCenter::const_iterator, Vff::Functional::construct_centers, 
        //!     Vff::functional::initialize_centers
        std::vector< t_Bond > bonds; 
        //! \brief Allow to relate i_atom_ pointer of Atomic_center::bonds to the correct
        //! periodic image with which the bond is made
        std::vector< math::rVector3d > translations;
        //! \brief A switch to say wether a bond is made directely or with a periodic
        //! image of an AtomicCenter
        std::vector< bool > do_translates;
        //! Crystal::Structure on which the Vff  functional is applied.
        Crystal :: TStructure<std::string> *structure;
        //! atomic index in Crystal::Structure::t_Atoms collection of AtomicCenter::structure
        types::t_unsigned index;

      public:
        mutable math::rVector3d gradient; //!< place holder to compute gradient

        //! \brief Default Constructor. 
        //! \param _str structure in which \a _e can be found
        //! \param _e atom to which this AtomicCenter relates
        //! \param _i index of _i in _str.atoms collection. Usefull for mpi processing
        AtomicCenter   (Crystal::TStructure<std::string> &_str, t_Atom &_e, types::t_unsigned _i) 
                     : i_atom_(&_e), structure(&_str), index(_i) {};
        //! \brief Copy Constructor
        //! \param[in] _c AtomicCenter object to copy
        AtomicCenter   ( const AtomicCenter &_c )
                      : i_atom_(_c.i_atom_), bonds(_c.bonds), translations(_c.translations), 
                        do_translates(_c.do_translates), structure(_c.structure),
                        gradient(0,0,0), index( _c.index) {} 

        //! \brief Adds _bond to  AtomicCenter::bonds if it is in first neighbor relationship
        //! \details This function returns -1 if \a _e is not a bond, and returns the
        //! number of bounds if it is. Note that poeriodic images of _bond are checked and 
        //! AtomicCenter::translations and AtomicCenter::do_translate are set accordingly.
        //! \param _bond AtomicCenter object for which first neighbor relationship is checked
        //! \param _cutoff distance below which first neighborness is implied
        types::t_int add_bond( t_BondRefd _bond, const types::t_real _cutoff  );

        //! Returns a AtomicCenter::const_iterator object pointing to the first bond
        const_iterator begin() const;
        //! Returns a AtomicCenter::const_iterator object pointing to the last bond
        const_iterator end() const;
        //! Returns the number of bonds
        types::t_unsigned size() const { return bonds.size(); }

        //! Sets the atomic position of the origin
        math::rVector3d& operator=(math::rVector3d& _vec) { i_atom_->pos = _vec; return _vec; }
        //! Translates the atomic position of the origin by _vec
        //! \param _vec Translation
        void operator+=(const math::rVector3d& _vec) { i_atom_->pos += _vec; }
        //! Translates the atomic position of the origin by -_vec
        //! \param _vec "Negative" Translation
        void operator-=(const math::rVector3d& _vec) { i_atom_->pos -= _vec; }
        //! Returns the atomic position of the origin
        operator math::rVector3d& () { return i_atom_->pos; }
        //! Returns the atomic position of the origin, constant format
        operator const math::rVector3d& () const { return i_atom_->pos; }
        //! Returns the atom at the origin
        t_Atom& atom() { return *i_atom_; }
        //! Returns the atom at the origin, constant format
        const t_Atom& atom() const { return *i_atom_; }
        //! Sets the gradient place-holder to the (0,0,0) vector
        void reset_gradient()
          { gradient[0] = 0; gradient[1] = 0; gradient[2] = 0; }
        //! Returns the index of this AtomicCenter in AtomicCenter::structure
        types::t_unsigned get_index() const { return index; }
    };

    //! \brief Iterator to travel along the bonds of an AtomicCenter
    //! \details Once a mesh of AtomicCenter objects is constructed, with the bonds
    //! forming the links in the mesh, an iterator is needed which will allow us to
    //! travel the mesh in any direction. AtomicCenter::const_iterator is built
    //! to iterates through the bonds of an AtomicCenter object. From there, one
    //! can easily travel throughout the mesh via first neighbor relationships.
    //! Dereferenced functions are performed by the end point of the bond.
    //! Underederenced functions are performed by the origin of the bond.
    class AtomicCenter :: const_iterator
    {
      //! The type of the atom  
      typedef Crystal::TStructure<std::string>::t_Atom  t_Atom;
      //! Type of pointer/iterator to the atomic center on the other side of the bond
      typedef AtomicCenter::t_Bond t_Bond;
      //! A reference to the of pointer/iterator to the atomic center on the other side of the bond
      typedef AtomicCenter::t_BondRefd t_BondRefd;
      //! \brief A constant reference to the of pointer/iterator to the atomic
      //         center on the other side of the bond
      typedef AtomicCenter::const_t_BondRefd const_t_BondRefd;

      protected:
        //! current origin of the bonds
        const AtomicCenter *parent;
        //! current bond being iterated
        std::vector< t_Bond > :: const_iterator i_bond;
        //! tranlation, if current bond is an atomic image
        std::vector< math::rVector3d > :: const_iterator i_translation;
        //! swtich for doing periodic imeage translation or not
        std::vector< bool > :: const_iterator i_do_translate;
      
      public:
        //! Constructor.
        const_iterator() {};
        //! \brief Constrcutor and Initialized
        //! \param _parent origin of the bond
        //! \param _is_begin
        //!                   - on true, initializes to first bond
        //!                   - on false, initializes to last bond
        //!                   .
        const_iterator   ( const AtomicCenter *_parent, bool _is_begin=true) 
                       : parent(_parent)
        {
          if ( _is_begin )
          {
            i_bond = parent->bonds.begin();
            i_translation = parent->translations.begin();
            i_do_translate = parent->do_translates.begin();
            return;
          }
          i_bond = parent->bonds.end();
#         ifdef LADA_DEBUG
            check();
#         endif
        }
        //! Copy Constrcutor 
        const_iterator   ( const const_iterator &_c ) 
                       : parent(_c.parent), i_bond(_c.i_bond), 
                         i_translation(_c.i_translation),
                         i_do_translate(_c.i_do_translate) {};
        //! Iterates through bonds
        void operator ++()
          { ++i_bond; ++i_translation; ++i_do_translate; }
        //! Iterates through bonds
        void operator --()
          { --i_bond; --i_translation; --i_do_translate; }
        //! Iterates through bonds
        void operator -=( types::t_int  n)
          { i_bond -= n; i_translation -= n; i_do_translate -= n; }
        //! Iterates through bonds
        void operator +=( types::t_int  n)
          { i_bond += n; i_translation += n; i_do_translate += n; }
        //! \brief Returns true if iterators are at same bond
        //! \param _i iterator agains which to check
        bool operator ==( const const_iterator &_i ) const
          { return _i.i_bond == i_bond; }
        //! \brief Returns true if iterators are \em not at same bond
        //! \param _i iterator agains which to check
        bool operator !=( const const_iterator &_i ) const
          { return _i.i_bond != i_bond; }
        //! \brief Returns the number of bonds separating this object  and _i in
        //!  bond array 
        //! \param _i iterator agains which to check
        types::t_int operator -( const const_iterator &_i )
          { return _i.i_bond - i_bond; }
        //! Dereferences the iterator: returns AtomicCenter to which bond points
        AtomicCenter& operator *()  { return *(*i_bond); }
        //! Dereferences the iterator: returns AtomicCenter to which bond points
        const AtomicCenter& operator *() const { return *(*i_bond); }
        //! Dereferences the iterator: returns AtomicCenter to which bond points
        const_t_BondRefd operator ->() const { return *i_bond; }
        //! Returns bond length squared
        types::t_real norm2() const
        {
#         ifdef LADA_DEBUG
            check_valid();
#         endif
          if ( not *i_do_translate )
            return (parent->i_atom_->pos - (*i_bond)->i_atom_->pos).squaredNorm();
          return (   parent->i_atom_->pos - (*i_bond)->i_atom_->pos
                   - parent->structure->cell * (*i_translation) ).squaredNorm();
        }
        //! \brief Returns bond vector
        math::rVector3d& vector( math::rVector3d &_hold )
        {
          _hold = (*i_bond)->i_atom_->pos - parent->i_atom_->pos ;
          if ( *i_do_translate )
            _hold += parent->structure->cell * (*i_translation);
          return _hold;
        }
        //! \brief Returns scalar product between bond vector and _b
        types::t_real scalar_product( const const_iterator &_b ) const
        {
          math::rVector3d a, b;
          if ( *i_do_translate )
            a =   (*i_bond)->i_atom_->pos - parent->i_atom_->pos 
                + parent->structure->cell * (*i_translation);
          else
            a = (*i_bond)->i_atom_->pos - parent->i_atom_->pos;
          if ( *_b.i_do_translate )
            b =   (*_b.i_bond)->i_atom_->pos - _b.parent->i_atom_->pos 
                + _b.parent->structure->cell * (*_b.i_translation);
          else
            b = (*_b.i_bond)->i_atom_->pos - _b.parent->i_atom_->pos;
          return a.dot(b);
        }
        //! Returns unique key for maps of bonds.
        std::string key() const
          { return bond_type(parent->atom().type, (*i_bond)->atom().type); }
        //! Returns unique key for maps of angles.
        std::string angle_key(AtomicCenter::const_iterator const &_other) const
        {
          return (_other.parent == parent) ? 
                 angle_type(parent->atom().type, (*i_bond)->atom().type, _other->atom().type):
                 angle_type(parent->atom().type, (*i_bond)->atom().type, _other.atom().type); 
        }
        //! \brief Returns the atom at the origin of range of bonds this iterator travels
        //! \see AtomicCenter::const_iterator::parent 
        t_Atom& atom() { return ((*i_bond)->atom()); }
        //! \brief Returns the atom at the origin of range of bonds this iterator travels
        //! \see AtomicCenter::const_iterator::parent 
        t_Atom const & atom() const { return ((*i_bond)->atom()); }
        //! \brief Translates a vector _v by periodic image of enpoint of bond
        //! \param _v vector to translate
        //! \param _cell unit-cell defining periodic image (can be different from
        //! AtomicCenter::structure by amount of minimized strain)
        void translate( math::rVector3d &_v, const math::rMatrix3d &_cell )
          { if( *i_do_translate ) _v += _cell * ( *i_translation ); }
        //! \brief Translates a vector _v by periodic image of enpoint of bond
        //! \param _v vector to translate
        void translate( math::rVector3d &_v )
          { if( *i_do_translate ) _v += parent->structure->cell * ( *i_translation ); }
#       ifdef LADA_DEBUG
          void check() const;
          void check_valid() const;
#       endif
    }; // end of const_iterator definition

    inline AtomicCenter::const_iterator AtomicCenter :: begin() const
      { return const_iterator( this ); }
    inline AtomicCenter::const_iterator AtomicCenter :: end() const
      { return const_iterator( this, false ); }

  } // namespace vff 
} // namespace LaDa

#endif // _VFF_FUNCTIONAL_H_

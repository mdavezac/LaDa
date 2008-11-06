//
//  Version: $Id$
//
#ifndef _VFF_ATOMIC_CENTER_H_
#define _VFF_ATOMIC_CENTER_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>

#include <tinyxml/tinyxml.h>

#include <crystal/atom.h>
#include <crystal/structure.h>
#include <crystal/lattice.h>
#include <opt/types.h>
#include <opt/function_base.h>
#include <opt/debug.h>

//! \brief Reimplements the Valence Force Field %Functional in c++
//! \details Vff, or Valence Force Field functional, is an empirical functional which
//! attempts to model the strain energy of a material from at most three-body
//! interactions. The there body interactions which are considered are
//! bond-stretching, change in bond-angles, and a combination of these two.
//! 
//! The implementation relies on a body-centered paradigm. In other words,
//! four classes have been created:
//!   - Vff::Functional is wrapper class and interface to Vff
//!   - Vff::Atomic_Center represent a single atom and lists its first neighbor relationships
//!   - Vff::Atomic_Center::const_iterator allows coders to travel
//!   through a collection of Vff::Atomic_Centera along first neighbor
//!   relationships.
//!   - Vff::Atomic_Functional computes the strain energy of one single atom, eg
//!   all three body terms in which a particular atom takes part.
//!   .
namespace Vff
{

  class Functional;

  //! \brief Represents a single Structure::t_Atom and its first neighbor relationships
  //! \details This class is meant to be used in conjunction with a list of
  //! Crystal::Structure::t_Atom, most likely in a Crystal::Structure. It contains a
  //! pointer, Atomic_Center::origin, which points a single Crystal::Structure::t_Atom. The
  //! first neighbor bonds of this atom are collected as vector of pointers to
  //! Vff::Atomic_Center objects in Atomic_Center::bonds. Since we are concerned
  //! with periodic structures, Atomic_Center::translations and
  //! Atomic_Center::bool record which periodic image of an atom
  //! Atomic_Center::bonds refer to.
  class Atomic_Center
  {
    friend class Functional;
    //! The type of the atom  
    typedef Crystal::Structure::t_Atom  t_Atom;
    //! The container of atomic centers. Defined here once and for all.
    typedef std::vector<Atomic_Center> t_Centers;
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
      t_Atom *origin; //!< The atom this object is addressing
      //! \brief Other Vff::Atomic_Center objects with which this one is in a bond-relationship
      //! \details Via bonds, a collection of Atomic_Center can be made into a tree,
      //! which can be travelled linearly, or through the first neighbor bonds,
      //! using Atomic_Center::const_iterator.
      //! \sa Atomic_Center::const_iterator, Vff::Functional::construct_centers, 
      //!     Vff::functional::initialize_centers
      std::vector< t_Bond > bonds; 
      //! \brief Allow to relate origin pointer of Atomic_center::bonds to the correct
      //! periodic image with which the bond is made
      std::vector< atat::rVector3d > translations;
      //! \brief A switch to say wether a bond is made directely or with a periodic
      //! image of an Atomic_Center
      std::vector< bool > do_translates;
      //! Crystal::Structure on which the Vff  functional is applied.
      Crystal :: Structure *structure;
      bool is_site_one; //!< helps determine the kind of atom this is
      bool is_site_one_two_species; //!< helps determine the kind of atom this is
      atat::rVector3d gradient; //!< place holder to compute gradient
      //! atomic index in Crystal::Structure::t_Atoms collection of Atomic_Center::structure
      types::t_unsigned index;

    public:
      //! \brief Default Constructor. 
      //! \param _str structure in which \a _e can be found
      //! \param _e atom to which this Atomic_Center relates
      //! \param _i index of _i in _str.atoms collection. Usefull for mpi processing
      Atomic_Center ( Crystal::Structure &_str, t_Atom &_e, types::t_unsigned _i);
      //! \brief Copy Constructor
      //! \param[in] _c Atomic_Center object to copy
      Atomic_Center   ( const Atomic_Center &_c )
                    : origin(_c.origin), bonds(_c.bonds), translations(_c.translations), 
                      do_translates(_c.do_translates), structure(_c.structure),
                      is_site_one(_c.is_site_one),
                      is_site_one_two_species( _c.is_site_one_two_species), 
                      gradient(0,0,0), index( _c.index) {} 

      //! \brief Returns the kind of atomic center this is
      //! \details In order to use the right coefficients in bond stretching and other
      //! interactions, we have to know the "kind" of
      //! functional this is. This will depend on the atomic specie of this
      //! atom and the encoding of Vff::Atomic_Functional array in
      //! Vff::Functional
      types::t_unsigned kind() const;

      //! \brief Adds _bond to  Atomic_Center::bonds if it is in first neighbor relationship
      //! \details This function returns -1 if \a _e is not a bond, and returns the
      //! number of bounds if it is. Note that poeriodic images of _bond are checked and 
      //! Atomic_Center::translations and Atomic_Center::do_translate are set accordingly.
      //! \param _bond Atomic_Center object for which first neighbor relationship is checked
      //! \param _cutoff distance below which first neighborness is implied
      types::t_int add_bond( t_BondRefd _bond, const types::t_real _cutoff  );

      //! Returns a Atomic_Center::const_iterator object pointing to the first bond
      const_iterator begin() const;
      //! Returns a Atomic_Center::const_iterator object pointing to the last bond
      const_iterator end() const;
      //! Returns the number of bonds
      types::t_unsigned size() const
        { return bonds.size(); }

      //! Sets the atomic position of the origin
      atat::rVector3d& operator=(atat::rVector3d& _vec)
        { origin->pos = _vec; return _vec; }
      //! Translates the atomic position of the origin by _vec
      //! \param _vec Translation
      void operator+=(const atat::rVector3d& _vec)
        { origin->pos += _vec; }
      //! Translates the atomic position of the origin by -_vec
      //! \param _vec "Negative" Translation
      void operator-=(const atat::rVector3d& _vec)
        { origin->pos -= _vec; }
      //! Returns the atomic position of the origin
      operator atat::rVector3d& ()
        { return origin->pos; }
      //! Returns the atomic position of the origin, constant format
      operator const atat::rVector3d& () const
        { return origin->pos; }
      //! Returns the atom at the origin
      t_Atom& Origin()
        { return *origin; }
      //! Returns the atom at the origin, constant format
      const t_Atom& Origin() const
        { return *origin; }
      //! Returns the gradient place holder
      atat::rVector3d& get_gradient()
        { return gradient; }
      //! Returns the gradient place holder, constant format
      const atat::rVector3d& get_gradient() const
        { return gradient; }
      //! Sets the gradient place-holder to the (0,0,0) vector
      void reset_gradient()
        { gradient[0] = 0; gradient[1] = 0; gradient[2] = 0; }
      //! Returns true if Atomic_Center is a site 1 in this lattice type
      bool site_one() const
        { return is_site_one; }
      //! Returns the index of this Atomic_Center in Atomic_Center::structure
      types::t_unsigned get_index() const
        { return index; }

    protected:
      //! \brief Returns the type of bond between this Atomic_Center and _bond
      //! \param _bond If this is not a bond, function will return result, not error!!
      //! \sa Atomic_Center::add_bond(), Atomic_Functional::add_bond(), Functional::Load()
      types::t_unsigned bond_kind( const Atomic_Center &_bond ) const;
  };

  //! \brief Iterator to travel along the bonds of an Atomic_Center
  //! \details Once a mesh of Atomic_Center objects is constructed, with the bonds
  //! forming the links in the mesh, an iterator is needed which will allow us to
  //! travel the mesh in any direction. Atomic_Center::const_iterator is built
  //! to iterates through the bonds of an Atomic_Center object. From there, one
  //! can easily travel throughout the mesh via first neighbor relationships.
  //! Dereferenced functions are performed by the end point of the bond.
  //! Underederenced functions are performed by the origin of the bond.
  class Atomic_Center :: const_iterator
  {
    //! The type of the atom  
    typedef Crystal::Structure::t_Atom  t_Atom;
    //! Type of pointer/iterator to the atomic center on the other side of the bond
    typedef Atomic_Center::t_Bond t_Bond;
    //! A reference to the of pointer/iterator to the atomic center on the other side of the bond
    typedef Atomic_Center::t_BondRefd t_BondRefd;
    //! \brief A constant reference to the of pointer/iterator to the atomic
    //         center on the other side of the bond
    typedef Atomic_Center::const_t_BondRefd const_t_BondRefd;

    protected:
      //! current origin of the bonds
      const Atomic_Center *parent;
      //! current bond being iterated
      std::vector< t_Bond > :: const_iterator i_bond;
      //! tranlation, if current bond is an atomic image
      std::vector< atat::rVector3d > :: const_iterator i_translation;
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
      const_iterator   ( const Atomic_Center *_parent, bool _is_begin=true) 
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
        __DOTRYDEBUGCODE( check();,
                             "error in constructor\n"
                          << " _parent != NULL " << (_parent!=NULL)
                          << " \n" << " _is_begin "
                          << _is_begin << "\n" )
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
      //! Dereferences the iterator: returns Atomic_Center to which bond points
      Atomic_Center& operator *()  { return *(*i_bond); }
      //! Dereferences the iterator: returns Atomic_Center to which bond points
      const Atomic_Center& operator *() const { return *(*i_bond); }
      //! Dereferences the iterator: returns Atomic_Center to which bond points
      const_t_BondRefd operator ->() const { return *i_bond; }
      //! Returns bond length squared
      types::t_real norm2() const
      {
        __DOTRYDEBUGCODE( check_valid();, "Invalid pointers in norm2()\n")
        if ( not *i_do_translate )
          return atat::norm2 ( parent->origin->pos - (*i_bond)->origin->pos );
        return atat::norm2( parent->origin->pos - (*i_bond)->origin->pos -
                            parent->structure->cell * (*i_translation) );
      }
      //! \brief Returns bond vector
      atat::rVector3d& vector( atat::rVector3d &_hold )
      {
        _hold = (*i_bond)->origin->pos - parent->origin->pos ;
        if ( *i_do_translate )
          _hold += parent->structure->cell * (*i_translation);
        return _hold;
      }
      //! \brief Returns scalar product between bond vector and _b
      types::t_real scalar_product( const const_iterator &_b ) const
      {
        atat::rVector3d a, b;
        if ( *i_do_translate )
          a =   (*i_bond)->origin->pos - parent->origin->pos 
              + parent->structure->cell * (*i_translation);
        else
          a = (*i_bond)->origin->pos - parent->origin->pos;
        if ( *_b.i_do_translate )
          b =   (*_b.i_bond)->origin->pos - _b.parent->origin->pos 
              + _b.parent->structure->cell * (*_b.i_translation);
        else
          b = (*_b.i_bond)->origin->pos - _b.parent->origin->pos;
        return a * b;
      }
      //! \brief Returns the kind of bond this is
      //! \see  Atomic_Center::bond_kind(), Atomic_Center::add_bond(),
      //!       Atomic_Functional::add_bond(), Functional::Load() 
      types::t_unsigned kind() const
        { return parent->bond_kind( *(*i_bond) ); }
      //! \brief Returns the atom at the origin of range of bonds this iterator travels
      //! \see Atomic_Center::const_iterator::parent 
      t_Atom& Origin()
        { return ((*i_bond)->Origin()); }
      //! \brief Translates a vector _v by periodic image of enpoint of bond
      //! \param _v vector to translate
      //! \param _cell unit-cell defining periodic image (can be different from
      //! Atomic_Center::structure by amount of minimized strain)
      void translate( atat::rVector3d &_v, const atat::rMatrix3d &_cell )
        { if( *i_do_translate ) _v += _cell * ( *i_translation ); }
      //! \brief Translates a vector _v by periodic image of enpoint of bond
      //! \param _v vector to translate
      void translate( atat::rVector3d &_v )
        { if( *i_do_translate ) _v += parent->structure->cell * ( *i_translation ); }
#ifdef _LADADEBUG
      void check() const;
      void check_valid() const;
#endif
  }; // end of const_iterator definition

  inline Atomic_Center::const_iterator Atomic_Center :: begin() const
    { return const_iterator( this ); }
  inline Atomic_Center::const_iterator Atomic_Center :: end() const
    { return const_iterator( this, false ); }

} // namespace vff 

#endif // _VFF_FUNCTIONAL_H_

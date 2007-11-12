//
//  Version: $Id$
//
#ifndef _LATTICE_H_
#define _LATTICE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <vector>
#include <fstream>
#include <string>
#include <complex>
#include <exception>


#include <tinyxml/tinyxml.h>
#include "opt/types.h"

#include "atat/findsym.h"
#include "atat/vectmac.h"
#include "atat/machdep.h"

#include "atom.h"


namespace Ising_CE {

  //! \brief Defines a lattice.
  //! \details A lattice is actually quite similar to a structure.
  //!          It contains the following elements:
  //!            - Lattice::cell defines the unit-cell in cartesian coordinate.
  //!            - Lattice::sites defines the each site in the unit-cell
  //!                             according to position in cartesian coordinate
  //!                             and their possible occupations.
  //!            - Lattice::scale defines the unit in which the cartesian coordinates are given.
  //!            - Lattice::space_group defines the the space-group of the lattice (How?
  //!            refer to <A
  //!            HREF="http://www.its.caltech.edu/~avdw/atat/manual/manual.html">
  //!            ATAT </A>)
  //!            .
  //!          Each site is an Atom_Type< std::vector<string> > and contains
  //!          both a cartesian position and set of possible occupations (in
  //!          strings of atomic symbol ).
  //
  //!          One of the main job of this class is to convert atoms with
  //!          atomic symbols into atoms with a numeric value. This only works
  //!          if each site can accept no more than two different %types of
  //!          atoms.
  //! \xmlinput A lattice can load itself (but not save to) from XML. See \ref TagLattice.
  //! \todo Make Structure into a template class? Then a lattice is no more
  //!       than a special structure...
  //! \warning Most indexing functions rely on Lattice::convert_type_index_to_real() and 
  //!          Lattice::convert_real_to_type_index(). In other words, they
  //!          won't work for sites which are not unaries or binaries...
  class Lattice
  {
    public:
      //! The type of each site.
      typedef Atom_Type< std::vector<std::string> > t_Site;
      //! The type of the collection of sites.
      typedef std::vector<t_Site> t_Sites;

    public:
      //! The unit-cell of the lattice in cartesian coordinates.
      atat::rMatrix3d cell;
      //! The collection of sites.
      t_Sites sites;
      //! The space-group of the lattice.
      atat::SpaceGroup space_group;
      //! The scale of the cartesian coordinates.
      types::t_real scale;

    public:
      //! Constructor.
      Lattice() {};
      //! Destructor.
      ~Lattice () {};

      //! Loads a lattice from XML.
      bool Load( const TiXmlElement &_element );

      //! \brief Computes the space-group of the lattice.
      //! \details Uses ATAT for this.
      void find_space_group();

      //! Returns the number of sites in the lattice.
      types::t_unsigned get_nb_sites() const
        { return sites.size(); } 
      //! Returns the number of possible occupations of site \a _i.
      types::t_unsigned get_nb_types( types::t_unsigned i ) const
        { return sites[i].type.size(); } 

      //! \brief Returns the site index of an atom at position \a _at.
      //! \details \a _at can be given modulo the unit-cell of the lattice.
      types::t_int get_atom_site_index( const atat::rVector3d &_at ) const;
      //! \brief Returns the site index of an atom with atomic symbol \a _at.
      //! \details Mores specifically, the first site in Lattice::sites which
      //!          accomodate \a _at is returned. This may not be the only site...
      types::t_int get_atom_site_index( const std::string &_at ) const;
      //! \brief Returns the index of the atomic occupation of \a _at.
      //! \details More specifically, returns the index of the string in member
      //!          Atom_Type< std::vector<std::string> >::type. The position
      //!          can be given modulo a lattice unit-cell.
      //! \warning works only as well as Lattice::convert_type_index_to_real().
      types::t_int get_atom_type_index( const Ising_CE::Atom &_at ) const;
      //! \brief Returns the index of the atomic occupation of \a _at.
      //! \details More specifically, returns the index of the string in member
      //!          Atom_Type< std::vector<std::string> >::type. 
      //!          If there are more than one site which can accomodate an \a
      //!          _at, the first the index in the first of the sites which can
      //!          is returned.
      types::t_int get_atom_type_index( const std::string &_at ) const;
      //! Returns the atomic symbol of \a _at.
      std::string get_atom_string( const Ising_CE::Atom &_at ) const;
      //! Returns the atomic symbol \a _i of site \a _s
      const std::string& get_atom_string( const unsigned _s, const unsigned _i ) const
        { return sites[_s].type[_i]; }
      //! Converts a atomic-symbol coded atom to an atom with a numeric value.
      bool convert_StrAtom_to_Atom( const Ising_CE::StrAtom &_in,
                                    Ising_CE::Atom &_out ) const;
      //! Converts a numerically coded atom to an atom with an atomic symbol.
      bool convert_Atom_to_StrAtom( const Ising_CE::Atom &_in,
                                    Ising_CE::StrAtom &_out ) const;
  
    protected:
      //! \brief Converts an index to a real value.
      //! \details Pretty much focused on binary sites.
      types::t_real convert_type_index_to_real( const types::t_unsigned _i ) const
        { return ( _i ) ? 1.0 : -1.0; }
      //! \brief Converts a real value to an index.
      //! \details Pretty much focused on binary and unary sites.
      types::t_unsigned convert_real_to_type_index( const types::t_real _r ) const
        { return ( std::abs(_r + 1.0) < atat::zero_tolerance ) ? 0 : 1; }
    public:
      //! \brief Converts the index \a _i of site \a _s to  a real value.
      //! \details Pretty much focused on binary and unary sites.
      types::t_real convert_type_index_to_real( const types::t_unsigned _s,
                                                const types::t_unsigned _i ) const
        { return ( sites[_s].type.size() == 1 ) ? -1: convert_type_index_to_real( _i ); }
      //! \brief Converts a real value to the index of site \a _s.
      //! \details Pretty much focused on binary and unary sites.
      types::t_unsigned convert_real_to_type_index( const types::t_unsigned _s,
                                                    const types::t_real _r ) const
        { return ( sites[_s].type.size() == 1 ) ? 0: convert_real_to_type_index( _r ); }
      //! Dumps the lattice to a stream. 
      void print_out (std::ostream &stream) const;
  };

  inline std::string Lattice::get_atom_string( const Ising_CE::Atom &_at ) const
  {
    types::t_int i = get_atom_site_index( _at );
    if ( i == -1 ) return "error";
    if ( get_nb_types(i) == 1 ) return sites[i].type[0];
    return ( std::abs( _at.type - 1.0 ) < atat::zero_tolerance ) ? 
        sites[i].type[0] : sites[i].type[1];  
  }

  //! Dumps a lattice to a stream.
  inline std::ostream& operator<<( std::ostream& _stream, const Ising_CE::Lattice& _lat )
    { _lat.print_out(_stream); return _stream; }
} // namespace Ising_CE
#endif

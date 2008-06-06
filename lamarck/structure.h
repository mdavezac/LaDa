//
//  Version: $Id$
//

#ifndef _ISING_CE_STRUCTURE_H_
#define _ISING_CE_STRUCTURE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include <ostream>
#include <fstream>
#include <string>
#include <complex>
#include <math.h>


#include <tinyxml/tinyxml.h>
#include <exception>

#include <opt/types.h>
#include <atat/vectmac.h>
#include <atat/xtalutil.h>
#include <atat/machdep.h>

#include "atom.h"
#include "lattice.h"

namespace Ising_CE {

  //! The atat structure type.
  typedef atat::Structure Atat_Structure;
  
  //! \brief Defines a periodic structure for a specific lattice.
  //! \details The definition is mostly sufficient and self contain. It include:
  //!          - Structure::cell the lattice cell
  //!          - Structure::atoms all atomic positions and occupations
  //!          - Structure::k_vecs all reciprocal-space vectors within the
  //!                                  first brillouin zone of the lattice
  //!                                  unit-cell.
  //!          - Structure::lattice a \a static pointer to the lattice. Static,
  //!                               because \e a \e priori we will work within
  //!                               a fixed lattice. Of course, this could
  //!                               change at some pointe, so it is well
  //!                               advised to use an instance to access the
  //!                               lattice.
  //!          - Structure::scale The scale in which the atomic positions and
  //!                             unit-cell are given. This should generally be
  //!                             the same as the scalar of the lattice.
  //!          - Structure::freeze which indicates which cartesian unit of the
  //!                              unit-cell should be frozen during
  //!                              optimization (say throug Vff).
  //!          .
  //! \xmlinput A structure can load and save itself to and from XML. See \ref
  //!           TagStructure. It is better to load the lattice first and the
  //!           structure(s) second.
  struct Structure
  {
    //! The atomic type
    typedef Atom_Type<types::t_real>  t_Atom;
    //! The type of the collection of atoms. 
    typedef std::vector< t_Atom >     t_Atoms;
    //! The reciprocal-space vector type
    typedef CAtom                     t_kAtom;
    //! The type of the collection of reciprocal-space vector.
    typedef std::vector< t_kAtom >    t_kAtoms;

    //! Freeze no coordinate of the unit-cell
    const static types::t_unsigned FREEZE_NONE =  0;
    //! Freeze the XX coordinate of the unit-cell ( Structure::cell(0,0) )
    const static types::t_unsigned FREEZE_XX   =  1;
    //! \brief Freeze the XY coordinate of the unit-cell ( Structure::cell(0,1) and
    //!        Structure::cell(1,0) )
    const static types::t_unsigned FREEZE_XY   =  2;
    //! \brief Freeze the XZ coordinate of the unit-cell ( Structure::cell(0,2) and
    //!        Structure::cell(2,0) )
    const static types::t_unsigned FREEZE_XZ   =  4;
    //! Freeze the YY coordinate of the unit-cell ( Structure::cell(2,2) )
    const static types::t_unsigned FREEZE_YY   =  8;
    //! \brief Freeze the YZ coordinate of the unit-cell ( Structure::cell(1,2) and
    //!        Structure::cell(2,1) )
    const static types::t_unsigned FREEZE_YZ   = 16;
    //! Freeze the ZZ coordinate of the unit-cell ( Structure::cell(2,2) )
    const static types::t_unsigned FREEZE_ZZ   = 32;
    //! Freeze all coordinates of the unit cell
    const static types::t_unsigned FREEZE_ALL  = 63;

    //! The unit-cell of the structure in cartesian coordinate.
    atat::rMatrix3d cell;
    //! The atomic position in cartesian unit and their occupation.
    std::vector< Atom_Type<types::t_real> > atoms;
    //! The reciprocal-space vector position in cartesian unit and their intensity.
    std::vector< CAtom > k_vecs;
    //! Just an old variable with the number of the structure in those NREL PI files.
    types::t_int Pi_name;
    //! The energy of the structure, whatever (and if) it is computed with.
    types::t_real energy;
    //! The frozen coordinates of the unit-cell.
    types::t_unsigned freeze;
    //! The scale in which cartesian units are given.
    types::t_real scale;
    //! A pointer to the lattice,
    static Lattice *lattice;

    public: 

    //! Constructor
    Structure() : Pi_name(0), energy(0), freeze(FREEZE_NONE) {};
    //! Constructor and Initializer
    Structure   ( const Atat_Structure &atat ) 
              : Pi_name(0), energy(0), freeze(FREEZE_NONE)
        { convert_from_ATAT( atat ); };
    //! Constructor. Loads itself from XML node \a _element.
    Structure   ( const TiXmlElement &_element )
              : Pi_name(0), energy(0), freeze(FREEZE_NONE)
      { Load( _element ); };
    //! Copy Constructor
    Structure   ( const Structure &_str )
              : cell(_str.cell), atoms(_str.atoms), k_vecs(_str.k_vecs),
                Pi_name(_str.Pi_name), energy(_str.energy), freeze(_str.freeze),
                scale( _str.scale ) {}
    //! Destructor.
    ~Structure () {};

    //! Converts an ATAT structure to an this class.
    void convert_from_ATAT (const Atat_Structure &atat);

    //! Converts an ATAT structure to an this class.
    void operator=( const Atat_Structure &atat )
     { convert_from_ATAT( atat ); };
    
    //! Prints a structure to a stream.
    void print_out( std::ostream &stream ) const;
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

    //! \brief Equates to \e geometic structure.
    //! \details The unit-cell and the atomic positions are compared, but not
    //!           the occupations.
    bool operator== (const Structure &_str ) const;
    
    //! \brief Copies the reciprocal-space vectors to a container.
    //! \details Since Atom_Type automatically returns reference to a
    //!          atat::rVector3d and its type, this routien can copy the full
    //!          reciprocal space vector, the positions only, or the
    //!          intensities only, depending on the type of
    //!          t_container::value_type.
    template<class t_container > void set_kvectors( const t_container &_container );

    //! Sets the site index of each atom according to Structure::lattice.
    bool set_site_indices();

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
    //! Serializes a structure.
    template<class ARCHIVE> void serialize(ARCHIVE & _ar, const unsigned int _version);
  };

  template< class ARCHIVE >
    void Structure :: serialize( ARCHIVE & _ar, const unsigned int _version)
    {
      _ar & cell;
      _ar & atoms;
      _ar & k_vecs;
      _ar & energy;
      _ar & freeze;
      _ar & scale;
    }

  //! \cond
  void  find_range( const atat::rMatrix3d &A, atat::iVector3d &kvec );
  //! \endcond
 
  //! Returns true if \a _a and \a _b are periodic equivalents of the unit-cell \a _cell.
  bool are_equivalent( const atat::rVector3d &_a,
                       const atat::rVector3d &_b,
                       const atat::rMatrix3d &_cell);
  //! Dumps a structure to a stream.
  inline std::ostream& operator<<( std::ostream& _stream, const Ising_CE::Structure& _struc )
    { _struc.print_out(_stream); return _stream; }

  //! compares two kvectors according to length and position.
  bool sort_kvec( const atat::rVector3d &_vec1, const atat::rVector3d &_vec2 );

  //! \brief Defines a %Fourier transform for a structure within a single-site lattice.
  //! \warning There is no reason this should work well in multiple-site lattices....
  struct Fourier
  {
    //! \brief From real to k space
    //! \details The first range is most likely some instance of
    //!          Ising_CE::Structure::t_Atoms. The second range is similarly an
    //!          instance of Ising_CE::Structure::t_kAtoms. The third range
    //! \param[in, out] _rfirst iterator to the first real space atom (of a type
    //!                similar to Ising_CE::Atom_Type)
    //! \param[in, out] _rend iterator to the last real space atom (of a type
    //!              similar to Ising_CE::Atom_Type)
    //! \param[in] _kfirst iterator to the first real space atom (of a type
    //!                similar to Ising_CE::Atom_Type< std::complex >)
    //! \param[in] _kend iterator to the last real space atom (of a type
    //!              similar to Ising_CE::Atom_Type< std::complex >)
    template<class T_R_IT, class T_K_IT>
    Fourier( T_R_IT _rfirst, T_R_IT _rend,
             T_K_IT _kfirst, T_K_IT _kend );
    //! \brief From k space to real space. 
    //! \details The first range is most likely some instance of
    //!          Ising_CE::Structure::t_Atoms. The second range is similarly an
    //!          instance of Ising_CE::Structure::t_kAtoms. The third range
    //!          should be iterators to std::complex.
    //! \pre The range [ \a _rout, \a _rout += \a _rfirst - \a _rend  ) should
    //!      be valid.
    //! \param[in] _rfirst iterator to the first real space atom (of a type
    //!                similar to Ising_CE::Atom_Type)
    //! \param[in] _rend iterator to the last real space atom (of a type
    //!              similar to Ising_CE::Atom_Type)
    //! \param[in] _kfirst iterator to the first real space atom (of a type
    //!                similar to Ising_CE::Atom_Type< std::complex >)
    //! \param[in] _kend iterator to the last real space atom (of a type
    //!              similar to Ising_CE::Atom_Type< std::complex >)
    //! \param[out] _rout iterator to the first complex real-space
    //!              occupation ( of std::complex type )
    template<class T_R_IT, class T_K_IT, class T_O_IT >
    Fourier( T_R_IT _rfirst, T_R_IT _rend,
             T_K_IT _kfirst, T_K_IT _kend,
             T_O_IT _rout ); // sets rvector values from kspace values
  };


  //! \cond
  template <class CONTAINER>
  void remove_equivalents( CONTAINER &_cont, const atat::rMatrix3d &_cell)
  {
    typename CONTAINER :: iterator i_vec = _cont.begin();
    typename CONTAINER :: iterator i_end = _cont.end();
    typename CONTAINER :: iterator i_which;

    while( i_vec != i_end )
    {
      i_which = i_vec+1;
      for ( ; i_which != i_end; i_which++ )
        if ( are_equivalent( *i_which, *i_vec, _cell ) )
          break;

      if ( i_which == i_end )
      { 
        ++i_vec;
        continue;
      }
      
      ( atat::norm2( (atat::rVector3d&) *i_vec ) < atat::norm2( (atat::rVector3d&) *i_which ) ) ? 
            _cont.erase(i_which): _cont.erase(i_vec);
      i_vec = _cont.begin();
      i_end = _cont.end();
    }

  }
  //! \endcond


  inline bool Structure :: operator== (const Structure &_str ) const
  {
    return     Fuzzy :: eq( cell(0,0), _str.cell(0,0) )
           and Fuzzy :: eq( cell(1,0), _str.cell(1,0) )
           and Fuzzy :: eq( cell(2,0), _str.cell(2,0) )
           and Fuzzy :: eq( cell(0,1), _str.cell(0,1) )
           and Fuzzy :: eq( cell(1,1), _str.cell(1,1) )
           and Fuzzy :: eq( cell(2,1), _str.cell(2,1) )
           and Fuzzy :: eq( cell(0,2), _str.cell(0,2) )
           and Fuzzy :: eq( cell(1,2), _str.cell(1,2) )
           and Fuzzy :: eq( cell(2,2), _str.cell(2,2) )
           and atoms == _str.atoms;
  }
    
  template<class t_container >
    void Structure :: set_kvectors( const t_container &_container )
    {
      typename t_container :: const_iterator i_kvec =  _container.begin();
      typename t_container :: const_iterator i_end =  _container.end();
      CAtom kvec;
      k_vecs.clear();
      k_vecs.reserve( _container.size() );
      for( ; i_kvec != i_end; ++i_kvec ) 
      {
        kvec = (*i_kvec);
        k_vecs.push_back( kvec );
      }
    }


  template<class T_R_IT, class T_K_IT>
  Fourier :: Fourier( T_R_IT _rfirst, T_R_IT _rend,
                      T_K_IT _kfirst, T_K_IT _kend )
  {
    const std::complex<types::t_real>
       imath(0, -2*3.1415926535897932384626433832795028841971693993751058208);
    
    for (; _kfirst != _kend; ++_kfirst)
    {
      _kfirst->type = std::complex<types::t_real>(0);
      for(T_R_IT i_r( _rfirst ); i_r != _rend; ++i_r )
      {
        _kfirst->type +=   exp( imath * ( i_r->pos[0] * _kfirst->pos[0] +
                                          i_r->pos[1] * _kfirst->pos[1] +
                                          i_r->pos[2] * _kfirst->pos[2] ) )
                         * i_r->type;
      }
    }
  }
  template<class T_R_IT, class T_K_IT, class T_O_IT >
  Fourier :: Fourier( T_R_IT _rfirst, T_R_IT _rend,
                      T_K_IT _kfirst, T_K_IT _kend,
                      T_O_IT _rout ) // sets rvector values from kspace values
  {
    const types::t_complex
       imath(0, 2*3.1415926535897932384626433832795028841971693993751058208);
    for (; _rfirst != _rend; ++_rfirst, ++_rout)
    {
      *_rout = types::t_complex(0,0);
      for(T_K_IT i_k=_kfirst; i_k != _kend; ++i_k)
      {
        *_rout +=   exp( imath * ( _rfirst->pos[0] * i_k->pos[0] +
                                   _rfirst->pos[1] * i_k->pos[1] +
                                   _rfirst->pos[2] * i_k->pos[2] ) )
                  * i_k->type;
      }
    }
  }




} // namespace Ising_CE


#endif

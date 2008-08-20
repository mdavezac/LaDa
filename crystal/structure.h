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
#include <boost/filesystem/path.hpp>


#include <tinyxml/tinyxml.h>
#include <exception>

#include <opt/types.h>
#include <atat/vectmac.h>
#include <atat/xtalutil.h>
#include <atat/machdep.h>

#include "atom.h"
#include "lattice.h"

namespace Crystal {

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
  template < class T_TYPE >
  struct TStructure
  {
    //! The atomic type
    typedef Atom_Type<T_TYPE>  t_Atom;
    //! The type of the collection of atoms. 
    typedef std::vector< t_Atom >     t_Atoms;

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
    std::vector< Atom_Type<T_TYPE> > atoms;
    //! Just an old variable with the number of the structure in those NREL PI files.
    std::string name;
    //! The energy of the structure, whatever (and if) it is computed with.
    types::t_real energy;
    //! Weight of structure in "some" set.
    types::t_real weight;
    //! The frozen coordinates of the unit-cell.
    types::t_unsigned freeze;
    //! The scale in which cartesian units are given.
    types::t_real scale;

    public: 

    //! Constructor
    TStructure() : name(""), energy(0), weight(1), freeze(FREEZE_NONE) {};
    //! Constructor. Loads itself from XML node \a _element.
    TStructure   ( const TiXmlElement &_element )
               : name(""), energy(0), weight(1), freeze(FREEZE_NONE)
      { Load( _element ); };
    //! Copy Constructor
    TStructure   ( const TStructure<T_TYPE> &_str )
               : cell(_str.cell), atoms(_str.atoms), name(_str.name), 
                 energy(_str.energy), weight(1), freeze(_str.freeze), scale( _str.scale ) {}
    //! Destructor.
    ~TStructure () {};

    //! Prints a structure to a stream.
    void print_out( std::ostream &stream ) const;
    //! Prints a structure to a string
    const std::string string() const
      { std::ostringstream sstr; print_out(sstr); return sstr.str(); }

    //! Loads a structure from XML.
    bool Load( const TiXmlElement &_element );
    //! Saves a structure to XML.
    void print_xml( TiXmlElement &_node ) const;

    //! \brief Equates to \e geometic structure.
    //! \details The unit-cell and the atomic positions are compared, but not
    //!           the occupations.
    bool operator== (const TStructure<T_TYPE> &_str ) const;
    
    //! Serializes a structure.
    template<class ARCHIVE> void serialize(ARCHIVE & _ar, const unsigned int _version);
  };
  
  struct Structure : public TStructure< types::t_real >
  {
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
    //! A pointer to the lattice,
    static Crystal::Lattice *lattice;

    public: 

    //! Constructor
    Structure() : TStructure<types::t_real>() {}
    //! Constructor and Initializer
    Structure   ( const Atat_Structure &atat ) 
              : TStructure<types::t_real>() { convert_from_ATAT( atat ); };
    //! Constructor. Loads itself from XML node \a _element.
    Structure   ( const TiXmlElement &_element )
              : TStructure<types::t_real>() { Load( _element ); };
    //! Copy Constructor
    Structure   ( const Structure &_str )
              : TStructure<types::t_real>( _str ), k_vecs(_str.k_vecs) {}
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

    //! \brief Copies the reciprocal-space vectors to a container.
    //! \details Since Atom_Type automatically returns reference to a
    //!          atat::rVector3d and its type, this routien can copy the full
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
    
    //! Sets the site index of each atom according to Structure::lattice.
    bool set_site_indices();

    //! Serializes a structure.
    template<class ARCHIVE> void serialize(ARCHIVE & _ar, const unsigned int _version);
  };

  //! Returns true if \a _a and \a _b are periodic equivalents of the unit-cell \a _cell.
  bool are_equivalent( const atat::rVector3d &_a,
                       const atat::rVector3d &_b,
                       const atat::rMatrix3d &_cell);
  //! Dumps a structure to a stream.
  inline std::ostream& operator<<( std::ostream& _stream, const Crystal::Structure& _struc )
    { _struc.print_out(_stream); return _stream; }

  //! compares two kvectors according to length and position.
  bool sort_kvec( const atat::rVector3d &_vec1, const atat::rVector3d &_vec2 );

  //! \brief Defines a %Fourier transform for a structure within a single-site lattice.
  //! \warning There is no reason this should work well in multiple-site lattices....
  struct Fourier
  {
    //! \brief From real to k space
    //! \details The first range is most likely some instance of
    //!          Crystal::Structure::t_Atoms. The second range is similarly an
    //!          instance of Crystal::Structure::t_kAtoms. The third range
    //! \param[in, out] _rfirst iterator to the first real space atom (of a type
    //!                similar to Crystal::Atom_Type)
    //! \param[in, out] _rend iterator to the last real space atom (of a type
    //!              similar to Crystal::Atom_Type)
    //! \param[in] _kfirst iterator to the first real space atom (of a type
    //!                similar to Crystal::Atom_Type< std::complex >)
    //! \param[in] _kend iterator to the last real space atom (of a type
    //!              similar to Crystal::Atom_Type< std::complex >)
    template<class T_R_IT, class T_K_IT>
    Fourier( T_R_IT _rfirst, T_R_IT _rend,
             T_K_IT _kfirst, T_K_IT _kend );
    //! \brief From k space to real space. 
    //! \details The first range is most likely some instance of
    //!          Crystal::Structure::t_Atoms. The second range is similarly an
    //!          instance of Crystal::Structure::t_kAtoms. The third range
    //!          should be iterators to std::complex.
    //! \pre The range [ \a _rout, \a _rout += \a _rfirst - \a _rend  ) should
    //!      be valid.
    //! \param[in] _rfirst iterator to the first real space atom (of a type
    //!                similar to Crystal::Atom_Type)
    //! \param[in] _rend iterator to the last real space atom (of a type
    //!              similar to Crystal::Atom_Type)
    //! \param[in] _kfirst iterator to the first real space atom (of a type
    //!                similar to Crystal::Atom_Type< std::complex >)
    //! \param[in] _kend iterator to the last real space atom (of a type
    //!              similar to Crystal::Atom_Type< std::complex >)
    //! \param[out] _rout iterator to the first complex real-space
    //!              occupation ( of std::complex type )
    template<class T_R_IT, class T_K_IT, class T_O_IT >
    Fourier( T_R_IT _rfirst, T_R_IT _rend,
             T_K_IT _kfirst, T_K_IT _kend,
             T_O_IT _rout ); // sets rvector values from kspace values
  };

  //! \brief fills in \a atoms member of an Crystal::Structure instance from the
  //!        cell-shape and the lattice.
  void fill_structure( Crystal::Structure &_str );
  //! Reads structure in NREL format.
  void read_structure( Structure &_struct, const boost::filesystem::path &_path );

  //! \brief reads lda energies and structures from NREL input files.
  //! \param[in] _path is the full or relative path to the "LDAs.dat" file.
  //!                  Structure files are expected to be found in the same
  //!                  directory as the LDAs.dat and is deduced from \a _path.
  //! \param[inout] _structures are added to this container.
  void read_ce_structures( const boost::filesystem::path &_dir,
                           std::vector<Crystal::Structure> &_structures );
  //! \brief Reads one structure from input file.
  //! \details Expects NREL pifile format. Returns true if succesfull, false if
  //!          reached eof.
  bool read_pifile_structure( std::istream &_sstr,
                              Crystal::Structure &_structure );

} // namespace Crystal

#include "structure.impl.h"


#endif

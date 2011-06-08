#ifndef _ISING_CE_STRUCTURE_H_
#define _ISING_CE_STRUCTURE_H_

#include "LaDaConfig.h"

#include <vector>
#include <ostream>
#include <fstream>
#include <string>
#include <complex>
#include <math.h>
#include <boost/filesystem/path.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/fusion/tuple/tuple_tie.hpp>


#include <tinyxml/tinyxml.h>
#include <exception>

#include <opt/types.h>
#include <opt/tinyxml.h>
#include <math/eigen.h>
#include <load_n_save/xpr/utilities.h>
#include <load_n_save/action/fusion.h>
#include <load_n_save/xpr/push_back.h>


#include "atom.h"
#include "lattice.h"


namespace LaDa 
{
  namespace Crystal
  {
    //! \cond
    class Lattice;
    //! \endcond

    namespace details
    {
      struct Structure
      {
        //! A pointer to the lattice,
        static Crystal::Lattice *lattice;
      };

    }

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
    struct TStructure : public details::Structure
    {
      friend class boost::serialization::access;
      //! The atomic type
      typedef Atom_Type<T_TYPE>  t_Atom;
      //! The type of the collection of atoms. 
      typedef std::vector< t_Atom >     t_Atoms;

      //! Tags to freeze cell coordinates.
      enum t_FreezeCell
      {
        FREEZE_NONE =   0, //!< Freeze no coordinate of the unit-cell
        FREEZE_XX   =   1, //!< Freeze (0,0)
        FREEZE_XY   =   2, //!< Freeze (0,1) 
        FREEZE_XZ   =   4, //!< Freeze (0,2)
        FREEZE_YX   =   8, //!< Freeze (0,1) 
        FREEZE_YY   =  16, //!< Freeze (1,1) 
        FREEZE_YZ   =  32, //!< Freeze (1,2) 
        FREEZE_ZX   =  64, //!< Freeze (2,2)
        FREEZE_ZY   = 128, //!< Freeze (2,2)
        FREEZE_ZZ   = 256, //!< Freeze (2,2)
        FREEZE_ALL  = 511, //!< Freeze all coordinates 
        FREEZE_A0   =  73, //!< Freeze all coordinates 
        FREEZE_A1   = 146, //!< Freeze all coordinates 
        FREEZE_A2   = 292, //!< Freeze all coordinates 
      };

      //! The unit-cell of the structure in cartesian coordinate.
      math::rMatrix3d cell;
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
                   energy(_str.energy), weight(_str.weight),
                   freeze(_str.freeze), scale( _str.scale ) {}
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
      
      //! Sets the site index of each atom according to Structure::lattice.
      bool set_site_indices();
      //! To load and save to xml-like input.
      template<class T_ARCHIVE>
        bool lns_access(T_ARCHIVE &_ar, load_n_save::version_type const _version) 
        {
          namespace lns = LaDa :: load_n_save;
          namespace bf = boost::fusion;
          std::map<std::string, LaDa::types::t_unsigned> freeze_map;
          freeze_map["none"] = FREEZE_NONE;
          freeze_map["a0"]   = FREEZE_A0;
          freeze_map["a1"]   = FREEZE_A1;
          freeze_map["a2"]   = FREEZE_A2;
          freeze_map["all"]  = FREEZE_ALL;
#         ifdef LADA_TIE
#            error LADA_TIE already defined.
#         endif
#         define LADA_TIE(i) bf::tie(cell(i,0), cell(i,1), cell(i,2))
#         ifdef LADA_TOE
#            error LADA_TOE already defined.
#         endif
#         define LADA_TOE(i) bf::tie(cell(0,i), cell(1,i), cell(2,i))
          lns::xpr::Section const seccell = lns::section("Cell") 
            << (
                    (    lns::option("r0", lns::tag=lns::required, lns::action=LADA_TIE(0))
                      && lns::option("r1", lns::tag=lns::required, lns::action=LADA_TIE(1))
                      && lns::option("r2", lns::tag=lns::required, lns::action=LADA_TIE(2)) )
                 || (    lns::option("a0", lns::tag=lns::required, lns::action=LADA_TOE(0))
                      && lns::option("a1", lns::tag=lns::required, lns::action=LADA_TOE(1))
                      && lns::option("a2", lns::tag=lns::required, lns::action=LADA_TOE(2))  )
               );
#         undef LADA_TIE
#         undef LADA_TOE
          t_Atom atom;
          lns::xpr::Section const section =
            lns::section("Structure")  
              << ( seccell  && lns::push_back(atoms) );
            // << lns::option("name", lns::action=name, lns::default_="")
            // << lns::option("energy", lns::action=energy, lns::default_=0)
            // << lns::option("weight", lns::action=weight, lns::default_=0)
            // << lns::option("freeze", lns::action=lns::enum_(freeze, freeze_map),
            //                lns::default_=FREEZE_NONE)
            // << lns::option("scale", lns::action=scale, lns::default_=1e0);
          return _ar & section;
        }

      private:
        //! Serializes a structure.
        template<class ARCHIVE> void serialize(ARCHIVE & _ar, const unsigned int _version);
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
      
//     //! Not found for some reason.
//     template<class T_ARCHIVE>
//       bool lns_access(const T_ARCHIVE &_ar, load_n_save::version_type const & _version) 
//         { return TStructure<types::t_real> :: lns_access(_ar, _version); }
      using TStructure<types::t_real> :: lns_access;

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

#ifndef LADA_CRYSTAL_LATTICE_H
#define LADA_CRYSTAL_LATTICE_H

#include "LaDaConfig.h"


#include <vector>
#include <fstream>
#include <string>
#include <complex>
#include <boost/shared_ptr.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/serialization/serialization.hpp>

#include <tinyxml/tinyxml.h>

#include <opt/debug.h>
#include <opt/types.h>
#include <opt/tinyxml.h>
#include <opt/path.h>

#include "atom.h"
#include "symmetry_operator.h"


namespace LaDa
{
  namespace Crystal 
  {
    //! \brief Refolds a periodic vector into the unit cell.
    //! \note Source code is in make_primitive.cc
    math::rVector3d into_cell( math::rVector3d const &_vec, 
                               math::rMatrix3d const &_cell, 
                               math::rMatrix3d const &_inv);
    
   
    //! \brief Refolds a periodic vector into the voronoi cell (eg first BZ or WignerSeitz).
    //! \note Source code is in make_primitive.cc
    math::rVector3d into_voronoi( math::rVector3d const &_vec, 
                                  math::rMatrix3d const &_cell, 
                                  math::rMatrix3d const &_inv);

    //! \brief Refolds a periodic vector into a cell centered around zero (in
    //!        fractional coordinates).
    //! \details Since the vector is refolded in fractional coordinates, it may
    //!          or may not be the vector with smallest norm. Use math::rVector3d
    //!          into_voronoi() to get the equivalent vector with smallest norm.
    //! \note Source code is in make_primitive.cc
    math::rVector3d zero_centered( math::rVector3d const &_vec, 
                                   math::rMatrix3d const &_cell, 
                                   math::rMatrix3d const &_inv);

    //! \brief Defines a lattice.
    //! \details A lattice is actually quite similar to a structure.
    //!          It contains the following elements:
    //!            - Lattice::cell defines the unit-cell in cartesian coordinate.
    //!            - Lattice::sites defines the each site in the unit-cell
    //!                             according to position in cartesian coordinate
    //!                             and their possible occupations.
    //!            - Lattice::scale defines the unit in which the cartesian coordinates are given.
    //!            - Lattice::space_group operations. Call Lattice::find_space()_group to initialize.
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
      friend class boost::serialization::access;
      public:
        //! The type of each site.
        typedef Atom_Type< std::vector<std::string> > t_Site;
        //! The type of the collection of sites.
        typedef std::vector<t_Site> t_Sites;

      public:
        //! Name of the lattice.
        std::string name;
        //! The unit-cell of the lattice in cartesian coordinates.
        math::rMatrix3d cell;
        //! The collection of sites.
        t_Sites sites;
        //! The space-group operations of the lattice.
        typedef std::vector<SymmetryOperator> t_SpaceGroup;
        //! The space-group operations of the lattice.
        t_SpaceGroup space_group;
        //! The scale of the cartesian coordinates.
        types::t_real scale;

        //! Constructor.
        Lattice() {};
        //! Copy Constructor.
        Lattice   ( Lattice const &_c ) 
                : name(_c.name), cell(_c.cell), sites(_c.sites),
                  space_group(_c.space_group), scale(_c.scale) {}
                     
        //! Destructor.
        ~Lattice () {};

        //! Loads a lattice from XML.
        bool Load( const TiXmlElement &_element );

        //! \brief Finds and stores space group operations.
        //! \param[in] _tol acceptable tolerance when determining symmetries.
        //!             -1 implies that types::tolerance is used.
        void find_space_group(types::t_real _tol = -1e0);

        //! Returns the number of sites in the lattice.
        types::t_unsigned get_nb_sites() const
          { return sites.size(); } 
        //! Returns the number of possible occupations of site \a _i.
        types::t_unsigned get_nb_types( types::t_unsigned i ) const
          { return sites[i].type.size(); } 

        //! \brief Returns the site index of an atom \a _at.
        //! \details \a _at can be given modulo the unit-cell of the lattice.
        template<class TTYPE>
        types::t_int get_atom_site_index( Atom_Type<TTYPE> &_at ) const
        {
          if( _at.site == -1 ) _at.site = get_atom_site_index( _at.pos ); 
          return _at.site;
        }
        //! \brief Returns the site index of an atom at position \a _at.
        //! \details \a _at can be given modulo the unit-cell of the lattice.
        types::t_int get_atom_site_index( const math::rVector3d &_at ) const;
        //! \brief Returns the site index of an atom with atomic symbol \a _at.
        //! \details Mores specifically, the first site in Lattice::sites which
        //!          accomodate \a _at is returned. This may not be the only site...
        types::t_int get_atom_site_index( const std::string &_at ) const;
        //! \brief Returns the index of the atomic occupation of \a _at.
        //! \details More specifically, returns the index of the string in member
        //!          Atom_Type< std::vector<std::string> >::type. The position
        //!          can be given modulo a lattice unit-cell.
        //! \warning works only as well as Lattice::convert_type_index_to_real().
        types::t_int get_atom_type_index( const Crystal::Atom &_at ) const;
        //! \brief Returns the index of the atomic occupation of \a _at.
        //! \details More specifically, returns the index of the string in member
        //!          Atom_Type< std::vector<std::string> >::type. 
        //!          If there are more than one site which can accomodate an \a
        //!          _at, the first the index in the first of the sites which can
        //!          is returned.
        types::t_int get_atom_type_index( const std::string &_at ) const;
        //! Returns the index of the atomic occupation of \a _at site \a _i.
        //! \details More specifically, returns the index of the string in member
        //!          Atom_Type< std::vector<std::string> >::type. 
        //!          If there are more than one site which can accomodate an \a
        //!          _at, the first the index in the first of the sites which can
        //!          is returned.
        types::t_int get_atom_type_index( const std::string &_at, types::t_unsigned _i ) const;
        //! Returns the atomic symbol of \a _at.
        std::string get_atom_string( const Crystal::Atom &_at ) const;
        //! Returns the atomic symbol \a _i of site \a _s
        const std::string& get_atom_string( const unsigned _s, const unsigned _i ) const
          { return sites[_s].type[_i]; }
        //! Converts a atomic-symbol coded atom to an atom with a numeric value.
        bool convert_StrAtom_to_Atom( const Crystal::StrAtom &_in,
                                      Crystal::Atom &_out ) const;
        //! Converts a numerically coded atom to an atom with an atomic symbol.
        bool convert_Atom_to_StrAtom( const Crystal::Atom &_in,
                                      Crystal::StrAtom &_out ) const;

        //! make lattice primitive. Returns true if it already is.
        bool make_primitive(types::t_real _tolerance = -1e0 );
    
      protected:
        //! \brief Converts an index to a real value.
        //! \details Pretty much focused on binary sites.
        types::t_real convert_type_index_to_real( const types::t_unsigned _i ) const
          { return ( _i ) ? 1.0 : -1.0; }
        //! \brief Converts a real value to an index.
        //! \details Pretty much focused on binary and unary sites.
        types::t_unsigned convert_real_to_type_index( const types::t_real _r ) const
          { return ( _r < 0e0 ) ? 0 : 1; }
      public:
        //! \brief Converts the index \a _i of site \a _s to  a real value.
        //! \details Pretty much focused on binary and unary sites.
        types::t_real convert_type_index_to_real( const types::t_unsigned _s,
                                                  const types::t_unsigned _i ) const
          { return ( sites[_s].type.size() == 1 ) ? 
                   -1: convert_type_index_to_real( _i ); }
        //! \brief Converts a real value to the index of site \a _s.
        //! \details Pretty much focused on binary and unary sites.
        types::t_unsigned convert_real_to_type_index( const types::t_unsigned _s,
                                                      const types::t_real _r ) const
          { return ( sites[_s].type.size() == 1 ) ?
                   0: convert_real_to_type_index( _r ); }
        //! Dumps the lattice to a stream. 
        void print_out (std::ostream &stream) const;
      private:
        //! Serializes a lattice.
        template<class ARCHIVE> void serialize(ARCHIVE & _ar, const unsigned int _version);
    };

    inline std::string Lattice::get_atom_string( const Crystal::Atom &_at ) const
    {
      types::t_int i;
#     ifdef LADA_DEBUG
        try
        {
#     endif
          i = get_atom_site_index( _at );
#     ifdef LADA_DEBUG
        }
        catch(...)
        {
          std::cerr << LADA_SPOT_ERROR
                    << "Caught error while converting string to numerical atom\n"
                    << _at << "\n";
          throw;
        }
#     endif
      if ( i == -1 ) return "error";
      if ( get_nb_types(i) == 1 ) return sites[i].type[0];
      return ( std::abs( _at.type - 1.0 ) < types::tolerance ) ? 
          sites[i].type[0] : sites[i].type[1];  
    }

    //! Dumps a lattice to a stream.
    inline std::ostream& operator<<( std::ostream& _stream, const Crystal::Lattice& _lat )
      { _lat.print_out(_stream); return _stream; }

    //! Load lattice with filename redirection.
    boost::shared_ptr< Crystal::Lattice > read_lattice( const TiXmlElement& _node );
    //! Load lattice from input file
    boost::shared_ptr< Crystal::Lattice > read_lattice( const boost::filesystem::path& _path );

    //! Reads lattice from input file \a _fpath in current directory, or in \a _dpath.
    boost::shared_ptr< Crystal::Lattice >
      read_lattice( const boost::filesystem::path &_fpath, 
                    const boost::filesystem::path &_dpath );
    //! Reads lattice from XML input in string format.
    boost::shared_ptr< Crystal::Lattice > read_lattice( const std::string& _string );
    //! Reads lattice from input file.
    boost::shared_ptr< Crystal::Lattice > read_lattice( const boost::filesystem::path& _path );

    //! Returns the number of species.
    types::t_unsigned nb_species( const Crystal::Lattice &_lattice );


    //! Returns true if a lattice has two sites and the two sites contain different species.
    bool lattice_has_same_species( const Lattice &_lattice );

    // This is inlined so that it gets recompiled with correct _CUBIC_CE_ or
    // _TETRAGONAL_CE_ each time.
    inline boost::shared_ptr< Crystal::Lattice > read_lattice( const std::string& _string )
    {
      LADA_TRY_BEGIN
      boost::shared_ptr< Crystal::Lattice > result( new Crystal::Lattice ); 

      TiXmlDocument doc;
      TiXmlDocument doc2;
      TiXmlHandle handle( &doc );
      doc.Parse( _string.c_str() );
      TiXmlElement *child = handle.FirstChild( "Job" )
                                  .FirstChild( "Lattice" ).Element();
      LADA_DO_NASSERT( not child, "Could not find Lattice in input.\n" )
      if( child->Attribute("filename") )
      {
        const boost::filesystem::path
          n( opt::expand_path( child->Attribute("filename") ) );
        LADA_DO_NASSERT( not boost::filesystem::exists( n ),
                    n.string() + " could not be found.\n" )
        std::string file;
        opt::read_xmlfile( n, doc2 );
        doc2.Parse( file.c_str() );
        LADA_DO_NASSERT( not doc2.FirstChild( "Job" ), "Job tag does not exist.\n" )
        LADA_DO_NASSERT( not doc2.FirstChild( "Job" )->FirstChildElement("Lattice"), 
                    "Lattice tag does not exist.\n" )
        child = doc2.FirstChild( "Job" )->FirstChildElement( "Lattice" );
        LADA_DO_NASSERT( not child, "Could not find Lattice in input.\n" )
      }
      LADA_DO_NASSERT( not result->Load(*child),
                  "Error while reading Lattice from input.\n")
#     if defined (_TETRAGONAL_CE_)
        // Only Constituent-Strain expects and space group determination
        // expect explicitely tetragonal lattice. 
        // Other expect a "cubic" lattice wich is implicitely tetragonal...
        // Historical bullshit from input structure files @ nrel.
        for( types::t_int i=0; i < 3; ++i ) 
          if( math::eq( result->cell(2,i), 0.5e0 ) )
            result->cell(2,i) = 0.6e0;
#     endif
      result->find_space_group();
#     if defined (_TETRAGONAL_CE_)
        // Only Constituent-Strain expects and space group determination
        // expect explicitely tetragonal lattice. 
        // Other expect a "cubic" lattice wich is implicitely tetragonal...
        // Historical bullshit from input structure files @ nrel.
        for( types::t_int i=0; i < 3; ++i ) 
          if( math::eq( result->cell(2,i), 0.6e0 ) )
            result->cell(2,i) = 0.5e0;
#     endif
      return result;
      LADA_TRY_END(, "Could not read lattice from input.\n" )
    }

    template< class ARCHIVE >
      void Lattice :: serialize( ARCHIVE & _ar, const unsigned int _version)
      {
        _ar & name;
        _ar & cell;
        _ar & sites;
        _ar & scale;
        _ar & space_group;
      }

  } // namespace Crystal
} // namespace LaDa
#endif

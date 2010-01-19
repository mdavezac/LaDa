//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<stdexcept> 
#include<set> 
#include<iterator> 
#include<algorithm> 

#include <Eigen/LU>

#include <print/manip.h>
#include <opt/tinyxml.h>
#include <math/misc.h>

#include "lattice.h"
#include "symmetry_operator.h"


namespace LaDa
{
  namespace Crystal {

         
    bool Lattice :: Load( const TiXmlElement &_element )
    {
      const TiXmlElement *child, *parent;
      double d; math::rVector3d vec;
      int i;

      std::string str = _element.Value();
      if ( str.compare("Lattice" ) != 0 )
        parent = _element.FirstChildElement("Lattice");
      else parent = &_element;
      __DOASSERT( not parent, "Could not find lattice tag on input.\n" )

      // reads cell first
      child = parent->FirstChildElement( "row" );
      for (i=0 ; child; ++i, child = child->NextSiblingElement("row") )
      {
        child->Attribute("x", &d); vec(0) = d;
        child->Attribute("y", &d); vec(1) = d;
        child->Attribute("z", &d); vec(2) = d;
        cell.row(i) = vec;
      }
      LADA_DOASSERT( i == 3, "Incorrect cell in lattice tag on input: " << i << *parent << "\n"; )

      // reads scale if it exist
      scale = 0;
      parent->Attribute("scale", &scale);
      if ( not scale ) scale = 1.0; 


      // then read sites
      child = parent->FirstChildElement("site" );
      __DOASSERT( not child, "Could not find site tag in lattice input.\n" )
      for (types::t_int i=0 ; child; i++, child = child->NextSiblingElement("site") )
      {
        t_Site site;
        child->Attribute("x", &d); site.pos(0) = d;
        child->Attribute("y", &d); site.pos(1) = d;
        child->Attribute("z", &d); site.pos(2) = d;
        const TiXmlElement *atom_xml = child->FirstChildElement("atom");
        site.type.clear();
        for ( ; atom_xml; atom_xml = atom_xml->NextSiblingElement("atom") )
        {
          std::string str = atom_xml->Attribute("type");
          site.type.push_back(str);
        }
        site.site = i;
        sites.push_back(site);
      }

      return true;
    }

    struct CompSites
    {
      CompSites( Lattice::t_Site::t_Type const &_types )
      {
        std::copy(_types.begin(), _types.end(), std::inserter(set_, set_.begin()));
      }
      CompSites( CompSites const &_c ): set_(_c.set_) {}
      bool operator()( Lattice::t_Site const &_site )
      {
        Lattice::t_Site::t_Type::const_iterator i_first = _site.type.begin();
        Lattice::t_Site::t_Type::const_iterator const i_end = _site.type.end();
        for(; i_first != i_end; ++i_first )
          if( set_.find(*i_first) == set_.end() ) break;
        return i_first == i_end;
      }

      std::set<Lattice::t_Site::t_Type::value_type> set_;
    };

    void Lattice :: find_space_group(types::t_real _tol)
    {
      boost::shared_ptr< t_SpaceGroup >
        syms(get_space_group_symmetries(*this, _tol));
      space_group = *syms;
    }

    types::t_int Lattice :: get_atom_site_index( const math::rVector3d &_at ) const
    {
      const math::rMatrix3d inv_cell = cell.inverse();
      std::vector< t_Site > :: const_iterator i_site = sites.begin(); 
      std::vector< t_Site > :: const_iterator i_end = sites.end(); 
      
      for (; i_site != i_end; ++i_site )
        if ( math::are_periodic_images(_at, i_site->pos, inv_cell) ) 
          return i_site - sites.begin();

      __THROW_ERROR("Could not find atomic site index!! " << _at << "\n" )
    }

    // returns  -1 on error
    types::t_int Lattice :: get_atom_site_index( const std::string &_at ) const
    {
      std::vector< t_Site > :: const_iterator i_site = sites.begin(); 
      std::vector< t_Site > :: const_iterator i_end = sites.end(); 
      std::vector< std::string > :: const_iterator i_type, i_type_end;
      for (; i_site != i_end; ++i_site )
      {
        i_type = i_site->type.begin();
        i_type_end = i_site->type.end();
        for (; i_type != i_type_end; ++i_type )
          if ( i_type->compare(_at) == 0 )
            return i_site - sites.begin();
      }

      __THROW_ERROR( "Could not find atomic site index!! " << _at << "\n" ) 
    }

    types::t_int Lattice :: get_atom_type_index( const Crystal :: Atom &_at ) const
    {
      const math::rMatrix3d inv_cell = cell.inverse();
      std::vector< t_Site > :: const_iterator i_site = sites.begin(); 
      std::vector< t_Site > :: const_iterator i_end = sites.end(); 
      
      for (; i_site != i_end; ++i_site )
        if (math::are_periodic_images(_at.pos, i_site->pos, inv_cell) ) 
        {
          if ( i_site->type.size() == 1 ) return 0;
          return convert_real_to_type_index( _at.type );
        }

      __THROW_ERROR("Could not find atomic site index!! " << _at << "\n" )
    }

   
    types::t_int Lattice :: get_atom_type_index( const std::string &_at ) const
    {
      std::vector< t_Site > :: const_iterator i_site = sites.begin(); 
      std::vector< t_Site > :: const_iterator i_end = sites.end(); 
      std::vector< std::string > :: const_iterator i_type, i_type_end;
      for (; i_site != i_end; ++i_site )
      {
        i_type = i_site->type.begin();
        i_type_end = i_site->type.end();
        for (; i_type != i_type_end; ++i_type )
          if ( i_type->compare(_at) == 0 )
            return i_type - i_site->type.begin();
      }

      __THROW_ERROR("Could not find atomic site index!! " << _at << "\n" )
    }
    
    types::t_int Lattice :: get_atom_type_index( const std::string &_at, types::t_unsigned _i ) const
    {
      __ASSERT( _i >= sites.size(), "Index out of range.\n" )
      const t_Site& site = sites[_i];

      const t_Site :: t_Type :: const_iterator i_first( site.type.begin() );
      const t_Site :: t_Type :: const_iterator i_last( site.type.end() );
      const t_Site :: t_Type :: const_iterator i_type
        = std::find( i_first, i_last, _at );
      __DOASSERT( i_type == i_last, "Could not find atomic specie " << _at << "\n" )

      return i_type - i_first;
    }
    
    bool Lattice :: convert_StrAtom_to_Atom( const Crystal::StrAtom &_in,
                                             Crystal::Atom &_out ) const
    {
      // on exit, _out is initialized, even is site and type not found
      _out.pos = _in.pos; 
      _out.freeze = _in.freeze;
      _out.type = 0;

      types::t_int site;
      __TRYDEBUGCODE( site = (_in.site > -1 ) ? _in.site : get_atom_site_index( _in.pos );,
                      "Caught error while converting numerical to string atom\n" )
      if( site < 0 ) return false;
      _out.site = site;
      if ( sites[site].type.size() == 1 )
        { _out.type = convert_type_index_to_real(0); return true; }
      __ASSERT( sites[site].type.size() > 2,
                   "More than two types at site " << sites[site].pos
                << "\nConversion should not work\n"; )

      types::t_int type;
      std::vector< std::string > :: const_iterator i_type, i_type_end;
      
      i_type = sites[site].type.begin();
      i_type_end = sites[site].type.end();
      for (type = 0; i_type != i_type_end; ++i_type, ++type )
        if ( i_type->compare(_in.type) == 0 ) break;
      if ( i_type == i_type_end ) // silent error, type not found
        return false;

      _out.type = convert_type_index_to_real( type );
      
      return true;
    }

    bool Lattice :: convert_Atom_to_StrAtom( const Crystal::Atom &_in,
                                             Crystal::StrAtom &_out ) const
    {
      // on exit, _out is initialized, even is site and type not found
      _out.pos = _in.pos; 
      _out.freeze = _in.freeze;
      _out.type = "";
      _out.site = _in.site;

      types::t_int site;
      __TRYDEBUGCODE( site = (_in.site > -1 ) ? _in.site : get_atom_site_index( _in.pos );,
                      "Caught error while converting numerical to string atom\n" )
      if( site < 0 ) return false;
      if ( sites[site].type.size() == 1 )
        { _out.type = sites[site].type[0]; return true; }
      __ASSERT( sites[site].type.size() > 2,
                   "More than two types at site " << sites[site].pos
                << "\nConversion should not work\n"; )

      types::t_unsigned type = convert_real_to_type_index( _in.type );
      _out.type = sites[site].type[type];
      
      return true;
    }

    void Lattice :: print_out (std::ostream &stream) const
    {
      stream << std::endl << " Lattice Unit-Cell\n"
             << cell << "\n"
             << " Lattice Sites " << "\n";
      t_Sites :: const_iterator i_site = sites.begin();
      t_Sites :: const_iterator i_end = sites.end();
      for( ; i_site != i_end; ++i_site )
      {
        stream << "  Position: ";
        i_site->print_out(stream);
        stream << "\n";
      }
    }

    bool lattice_has_same_species( const Lattice &_lattice )
    {
      if( _lattice.sites.size() != 2 ) return false;
      const Lattice :: t_Site :: t_Type &site1 = _lattice.sites[0].type;
      const Lattice :: t_Site :: t_Type &site2 = _lattice.sites[1].type;

      foreach( const Lattice :: t_Site :: t_Type :: value_type& type, site1 )
        if( std::find( site2.begin(), site2.end(), type ) == site2.end() ) return false;

      return true;
    }

    types::t_unsigned nb_species( const Lattice &_lattice )
    {
      std::set< std::string > set;
      foreach( const Lattice::t_Site site, _lattice.sites )
        std::copy( site.type.begin(), site.type.end(), std::inserter( set, set.begin() ) );
      return set.size();
    }

    boost::shared_ptr< Crystal::Lattice >
      read_lattice( const boost::filesystem::path &_fpath, 
                    const boost::filesystem::path &_dpath )
      { 
        boost::filesystem::path filepath( _fpath );
        if( not boost::filesystem::exists( filepath ) )
          filepath = _dpath / _fpath;
        __DOASSERT( not boost::filesystem::exists( filepath ), 
                       "Could not find "<< _fpath 
                    << " in current directory, nor in " <<  _dpath )
        return read_lattice( filepath );
      }
    boost::shared_ptr< Crystal::Lattice > read_lattice( const boost::filesystem::path& _path )
    {
      __TRYBEGIN
        std::string input_filestream;
        opt::read_xmlfile( _path, input_filestream );
        return read_lattice( input_filestream );
      __TRYEND(, "Could not read lattice from input.\n" )
    }

    boost::shared_ptr< Crystal::Lattice > read_lattice( const TiXmlElement& _node )
    {
      __TRYBEGIN
        boost::shared_ptr< Crystal::Lattice > lattice( new Lattice );
        opt::read_tag( *lattice, _node, "Lattice" );
        return lattice;
      __TRYEND(,"Could not real lattice from input.\n" )
    }


    //  bool make_primitive() is to be found in make_primitive.cc

  } // namespace Crystal

} // namespace LaDa

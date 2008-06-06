//
//  Version: $Id$
//
#include<stdexcept> 

#include <atat/array.h>
#include <atat/misc.h>
#include <opt/ndim_iterator.h>

#include "lattice.h"


namespace Ising_CE {

       
  bool Lattice :: Load( const TiXmlElement &_element )
  {
    const TiXmlElement *child, *parent;
    double d; atat::rVector3d vec;
    int i;

    std::string str = _element.Value();
    if ( str.compare("Lattice" ) != 0 )
      parent = _element.FirstChildElement("Lattice");
    else
      parent = &_element;
    __DOASSERT( not parent, "Could not find lattice tag on input.\n" )

    // reads cell first
    child = parent->FirstChildElement( "column" );
    for (i=0 ; child; i++, child = child->NextSiblingElement("column") )
    {
      child->Attribute("x", &d); vec(0) = d;
      child->Attribute("y", &d); vec(1) = d;
      child->Attribute("z", &d); vec(2) = d;
      cell.set_column(i, vec);
    }
    __DOASSERT( i != 3, "Incorrect cell in lattice tag on input.\n" )

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

  // finds space group and symmetry operations using atat library
  // needs some hooking in
  // also wraps inside cell
  void Lattice :: find_space_group()
  {
    space_group.cell=cell;
    atat::Array<atat::rVector3d> atom_pos( sites.size() );
    atat::Array<types::t_int> atom_type( sites.size() );
    std::vector<t_Site> :: iterator i_site = sites.begin();
    std::vector<t_Site> :: iterator i_end = sites.end();

    // copies to an atat::Array
    for(types::t_unsigned n=0; i_site != i_end; ++i_site, ++n)
    {
      atom_pos[n] = i_site->pos;
      atom_type[n] = n;
    }
    atat::wrap_inside_cell( &atom_pos, atom_pos, cell );
    // recopies wrapped atoms back to this object
    i_site = sites.begin();
    for(types::t_unsigned n=0; i_site != i_end; ++i_site, ++n)
      i_site->pos = atom_pos[n];

    // finally, looks for space group
    find_spacegroup(&space_group.point_op, &space_group.trans,
                    cell, atom_pos, atom_type);
    if (contains_pure_translations(space_group.point_op,space_group.trans)) 
      std::cerr << "Warning: unit cell is not primitive." << std::endl;
    // Makes sure that translations are not by some integer combination of the
    // unit cell.
    for( types::t_int i = 0; i < space_group.trans.getSize(); ++i )
    {
      atat::rVector3d &trans = space_group.trans(i);
      if( Fuzzy::eq( atat::norm2( trans ), types::t_real(0) ) ) continue;
      atat::rVector3d zeroed = (!cell) * trans;
      zeroed[0] = zeroed[0] - std::floor( zeroed[0] + 0.1 ); 
      zeroed[1] = zeroed[1] - std::floor( zeroed[1] + 0.1 ); 
      zeroed[2] = zeroed[2] - std::floor( zeroed[2] + 0.1 ); 
      if( Fuzzy::eq( zeroed[0], types::t_real(0) ) ) zeroed[0] = types::t_real( 0 );
      if( Fuzzy::eq( zeroed[1], types::t_real(0) ) ) zeroed[1] = types::t_real( 0 );
      if( Fuzzy::eq( zeroed[2], types::t_real(0) ) ) zeroed[2] = types::t_real( 0 );
      trans = zeroed;
    }
    
  }

  types::t_int Lattice :: get_atom_site_index( const atat::rVector3d &_at ) const
  {
    const atat::rMatrix3d inv_cell = !cell;
    std::vector< t_Site > :: const_iterator i_site = sites.begin(); 
    std::vector< t_Site > :: const_iterator i_end = sites.end(); 
    
    for (; i_site != i_end; ++i_site )
      if ( atat::equivalent_mod_cell(_at, i_site->pos,inv_cell) ) 
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

  types::t_int Lattice :: get_atom_type_index( const Ising_CE :: Atom &_at ) const
  {
    const atat::rMatrix3d inv_cell = !cell;
    std::vector< t_Site > :: const_iterator i_site = sites.begin(); 
    std::vector< t_Site > :: const_iterator i_end = sites.end(); 
    
    for (; i_site != i_end; ++i_site )
      if ( atat::equivalent_mod_cell(_at.pos, i_site->pos,inv_cell) ) 
      {
        if ( i_site->type.size() == 1 ) return 0;
        return convert_real_to_type_index( _at.type );
      }

    __THROW_ERROR("Could not find atomic site index!! " << _at << "\n" )
  }

  // returns  -1 on error
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
  
  bool Lattice :: convert_StrAtom_to_Atom( const Ising_CE::StrAtom &_in,
                                           Ising_CE::Atom &_out ) const
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

  bool Lattice :: convert_Atom_to_StrAtom( const Ising_CE::Atom &_in,
                                           Ising_CE::StrAtom &_out ) const
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
    stream << std::endl << " Lattice Unit-Cell " << std::endl;
    stream << cell;
    
    stream << " Lattice Sites " << std::endl;
    t_Sites :: const_iterator i_site = sites.begin();
    t_Sites :: const_iterator i_end = sites.end();
    for( ; i_site != i_end; ++i_site )
    {
      stream << "  Position: ";
      i_site->print_out(stream);
      stream << "\n";
    }
  }

  // refold by one vector
  void refold( atat::rVector3d &vec, const atat::rMatrix3d &lat )
  {
    opt::Ndim_Iterator<types::t_int, std::less_equal<types::t_int> > i_cell;
    atat::rVector3d hold = vec;
    atat::rVector3d compute;
    atat::rVector3d current = vec;
    types::t_real norm_c = norm2(vec);

    i_cell.add(-2,2);
    i_cell.add(-2,2);
    i_cell.add(-2,2);

    do
    {
      compute(0) = (types::t_real) i_cell.access(0);
      compute(1) = (types::t_real) i_cell.access(1);
      compute(2) = (types::t_real) i_cell.access(2);

      vec = hold + lat*compute;
      if ( norm2( vec ) < norm_c ) 
      {
        current = vec;
        norm_c = norm2(vec);
      }

    } while ( (++i_cell) );

    vec = current;
  }

} // namespace Ising_CE


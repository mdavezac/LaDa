#include "atat/array.h"
#include "atat/misc.h"
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
    if ( not parent )
      return false;

    // reads cell first
    child = parent->FirstChildElement( "column" );
    for (i=0 ; child; i++, child = child->NextSiblingElement("column") )
    {
      child->Attribute("x", &d); vec(0) = d;
      child->Attribute("y", &d); vec(1) = d;
      child->Attribute("z", &d); vec(2) = d;
      cell.set_column(i, vec);
    }
    if ( i != 3 )
      return false;

    // reads scale if it exist
    scale = 0;
    parent->Attribute("scale", &scale);
    if ( not scale ) scale = 1.0; 


    // then read sites
    child = parent->FirstChildElement("site" );
    if ( not child )
      return false;
    for (types::t_int i=0 ; child; i++, child = child->NextSiblingElement("site") )
    {
      t_Site site;
      child->Attribute("x", &d); site.pos(0) = d;
      child->Attribute("y", &d); site.pos(1) = d;
      child->Attribute("z", &d); site.pos(2) = d;
      const TiXmlElement *atom_xml = child->FirstChildElement("atom");
      site.type.clear();
      for ( ; atom_xml; i++, atom_xml = atom_xml->NextSiblingElement("atom") )
      {
        std::string str = atom_xml->Attribute("type");
        site.type.push_back(str);
      }
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
  }

  types::t_int Lattice :: get_atom_site_index( const atat::rVector3d &_at ) const
  {
    const atat::rMatrix3d inv_cell = !cell;
    std::vector< t_Site > :: const_iterator i_site = sites.begin(); 
    std::vector< t_Site > :: const_iterator i_end = sites.end(); 
    
    for (; i_site != i_end; ++i_site )
      if ( atat::equivalent_mod_cell(_at, i_site->pos,inv_cell) ) 
        return i_site - sites.begin();

    std::cerr << "Could not find equivalent site!! " << std::endl;
    return -1;
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

    std::cerr << "Could not find atom site " << _at << "!!" << std::endl;
    return -1;
  }

  types::t_int Lattice :: get_atom_type_index( const Ising_CE :: Atom &_at ) const
  {
    const atat::rMatrix3d inv_cell = !cell;
    std::vector< t_Site > :: const_iterator i_site = sites.begin(); 
    std::vector< t_Site > :: const_iterator i_end = sites.end(); 
    
    for (; i_site != i_end; ++i_site )
      if ( atat::equivalent_mod_cell(_at.pos, i_site->pos,inv_cell) ) 
      {
        if ( i_site->type.size() == 1 )
          return 0;
        return convert_real_to_type_index( _at.type );
      }

    std::cerr << "Could not find equivalent site!! " << std::endl;
    return -1;
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
    std::cerr << "Could not find atom site " << _at << "!!" << std::endl;
    return -1;
  }
  
  bool Lattice :: convert_StrAtom_to_Atom( const Ising_CE::StrAtom &_in,
                                           Ising_CE::Atom &_out ) const
  {
    // on exit, _out is initialized, even is site and type not found
    _out.pos = _in.pos; 
    _out.freeze = _in.freeze;
    _out.type = 0;
    _out.site = _in.site;

    types::t_int site = (_in.site > -1 ) ? _in.site : get_atom_site_index( _in.pos );
    if ( sites[site].type.size() == 1 )
      { _out.type = convert_type_index_to_real(0); return true; }
    if ( sites[site].type.size() > 2 )
      std::cerr << "More than three types at site " << sites[site].pos
                << std::endl << "Conversion should not work " << std::endl;

    types::t_int type;
    std::vector< std::string > :: const_iterator i_type, i_type_end;
    
    i_type = sites[site].type.begin();
    i_type_end = sites[site].type.end();
    for (type = 0; i_type != i_type_end; ++i_type, ++type )
      if ( i_type->compare(_in.type) == 0 )
        break;
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

    types::t_int site = (_in.site > -1 ) ? _in.site : get_atom_site_index( _in.pos );
    if ( sites[site].type.size() == 1 )
      { _out.type = sites[site].type[0]; return true; }
    if ( sites[site].type.size() > 2 )
      std::cerr << "More than three types at site " << sites[site].pos
                << std::endl << "Conversion should not work " << std::endl;

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
      stream << "  Position: " << i_site->pos[0] << " "
             << i_site->pos[1] << " "
             << i_site->pos[2] << " ";
      for( types::t_unsigned i = 0; i < i_site->type.size(); ++i )
        stream << i_site->type[i] << " ";
      stream << std::endl;
    }
  }

} // namespace Ising_CE

#ifdef _MPI
#include<mpi/mpi_object.h>

namespace mpi
{
  template<>
  bool BroadCast :: serialize<atat::SpaceGroup>( atat::SpaceGroup &_sg )
  {
    // first copies cell vectors
    if ( not serialize( _sg.cell.x[0], (_sg.cell.x[0]+3) ) ) return false;
    if ( not serialize( _sg.cell.x[1], (_sg.cell.x[1]+3) ) ) return false;
    if ( not serialize( _sg.cell.x[2], (_sg.cell.x[2]+3) ) ) return false;

    // then serializes point operations
    types::t_int n = _sg.point_op.get_size();
    if( not serialize( n ) ) return false;
    if ( stage == COPYING_FROM_HERE )
      _sg.point_op.resize(n);
    for( types::t_int i = 0; i < n ; ++i )
    {
      if ( not serialize( _sg.point_op[i].x[0], (_sg.point_op[i].x[0]+3) ) ) return false;
      if ( not serialize( _sg.point_op[i].x[1], (_sg.point_op[i].x[1]+3) ) ) return false;
      if ( not serialize( _sg.point_op[i].x[2], (_sg.point_op[i].x[2]+3) ) ) return false;
    }
    
    // finally serializes translations operations
    n = _sg.trans.get_size();
    if( not serialize( n ) ) return false;
    if ( stage == COPYING_FROM_HERE )
      _sg.trans.resize(n);
    for( types::t_int i = 0; i < n ; ++i )
      if ( not serialize( _sg.trans[i].x, (_sg.trans[i].x+3) ) ) return false;

    return true;
  }

  template<>
  bool BroadCast :: serialize<Ising_CE::Lattice>( Ising_CE::Lattice &_latt )
  {
    // first copies cell vectors
    if ( not serialize( _latt.cell.x[0], (_latt.cell.x[0]+3) ) ) return false;
    if ( not serialize( _latt.cell.x[1], (_latt.cell.x[1]+3) ) ) return false;
    if ( not serialize( _latt.cell.x[2], (_latt.cell.x[2]+3) ) ) return false;
    
    // copies scale
    if ( not serialize( _latt.scale ) ) return false;

    // copies site size
    types::t_int n = _latt.sites.size();
    if( not serialize( n ) ) return false;
    if ( stage == COPYING_FROM_HERE )
      _latt.sites.resize(n);

    // copies sites
    Ising_CE::Lattice::t_Sites :: iterator i_site = _latt.sites.begin();
    Ising_CE::Lattice::t_Sites :: iterator i_site_end = _latt.sites.end();
    for(; i_site != i_site_end; ++i_site )
    {
      if ( not serialize( i_site->pos.x, i_site->pos.x+3 ) ) return false;
      if ( not serialize( i_site->freeze ) ) return false;
      n = i_site->type.size();
      if( not serialize( n ) ) return false;
      if ( stage == COPYING_FROM_HERE )
        i_site->type.resize(n);
      std::vector<std::string> :: iterator i_str = i_site->type.begin();
      std::vector<std::string> :: iterator i_str_end = i_site->type.end();
      for(; i_str != i_str_end; ++i_str )
        if( not serialize( *i_str ) ) return false;
    }

    // finally serializes spacegroup
    if (not serialize( _latt.space_group ) ) return false;

    return true;
  }
}

#endif

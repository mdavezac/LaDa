//
//  Version: $Id$
//
#include <algorithm>

#include <opt/ndim_iterator.h>
#include <opt/traits.h>

#include <atat/misc.h>
#include <physics/physics.h>

#include "structure.h"




namespace Crystal {

  bool sort_kvec( const atat::rVector3d &_vec1, const atat::rVector3d &_vec2 )
  {
    types::t_real a = atat::norm2( _vec1 );
    types::t_real b = atat::norm2( _vec2 );
    if ( std::abs( a - b ) > types::tolerance ) return a < b;
    a = _vec1[0]; b = _vec2[0];
    if ( std::abs( a - b ) > types::tolerance ) return a < b;
    a = _vec1[1]; b = _vec2[1];
    if ( std::abs( a - b ) > types::tolerance ) return a < b;
    a = _vec1[2]; b = _vec2[2];
    if ( std::abs( a - b ) > types::tolerance ) return a < b;
    return false;
  }

  using atat::rndseed;
  Crystal::Lattice* Structure :: lattice = NULL;

  void Structure :: convert_from_ATAT ( const Atat_Structure &atat  )
  {
    cell = atat.cell;
    for ( types::t_int i = 0; i < atat.atom_pos.get_size(); i++ )
      atoms.push_back( Atom_Type<types::t_real>( atat.atom_pos[i],
                             static_cast<types::t_real>(atat.atom_type[i]) ) );
  }
                                                                 
  void Structure :: print_out (std::ostream &stream) const
  {
    stream << "\n Structure, scale: " << scale << ", Volume: "
           << atat::det( cell )
           << ", Cell\n"
           << std::fixed << std::setprecision(5)
           << "   " << std::setw(9) << cell(0,0)
           << "   " << std::setw(9) << cell(1,0)
           << "   " << std::setw(9) << cell(2,0) << "\n"
           << "   " << std::setw(9) << cell(0,1)
           << "   " << std::setw(9) << cell(1,1)
           << "   " << std::setw(9) << cell(2,1) << "\n"
           << "   " << std::setw(9) << cell(0,2)
           << "   " << std::setw(9) << cell(1,2)
           << "   " << std::setw(9) << cell(2,2) << "\n"
           << "\n   Atoms:\n";
            
    
    stream << " Structure atoms " << std::endl;
    std::vector<Atom> :: const_iterator i_atom = atoms.begin();
    std::vector<Atom> :: const_iterator i_end = atoms.end();
    for( ; i_atom != i_end; ++i_atom )
    {
      stream << "  Position: ";
      StrAtom stratom;
      if( lattice and lattice->convert_Atom_to_StrAtom( *i_atom, stratom ) )
        stratom.print_out(stream); 
      else i_atom->print_out(stream); 
      stream << std::endl;
    }
    if ( not k_vecs.size() ) return;
    stream << " Structure K vectors " << std::endl;
    t_kAtoms :: const_iterator i_kvec = k_vecs.begin();
    t_kAtoms :: const_iterator i_kvec_end = k_vecs.end();
    for( ; i_kvec != i_kvec_end; ++i_kvec )
      stream << "  Kvec: " 
             << std::fixed << std::setprecision(5)
             << *i_kvec << "\n";
  }


  void Structure :: set_atom_types( const std::vector<types::t_real> &types)
  {
    __ASSERT( types.size() != atoms.size(), "Inequivalent vectors.\n" )

    std::vector<types::t_real> :: const_iterator i_type = types.begin();
    std::vector<types::t_real> :: const_iterator i_type_last = types.end();
    std::vector<Atom> :: iterator i_atom = atoms.begin();
    std::vector<Atom> :: iterator i_atom_last = atoms.end();

    for ( ; i_type != i_type_last and i_atom != i_atom_last;
          ++i_type, ++i_atom)
        i_atom->type = *i_type;
  }

  void Structure :: get_atom_types( std::vector<types::t_real> &types) const
  {
    __ASSERT( types.size() != atoms.size(), "Inequivalent vectors.\n" )

    std::vector<types::t_real> :: iterator i_type = types.begin();
    std::vector<types::t_real> :: iterator i_type_last = types.end();
    std::vector<Atom> :: const_iterator i_atom = atoms.begin();
    std::vector<Atom> :: const_iterator i_atom_last = atoms.end();

    for ( ; i_type != i_type_last and i_atom != i_atom_last;
          ++i_type, ++i_atom)
      *i_type = i_atom->type;
  }

  types::t_real Structure :: get_concentration() const
  {
    std::vector<Atom> :: const_iterator i_atom = atoms.begin();
    std::vector<Atom> :: const_iterator i_atom_last = atoms.end();
    types::t_real result = 0;
    for ( ; i_atom != i_atom_last; i_atom++)
      result += i_atom->type;

    return ( result / ( (types::t_real) atoms.size() ) );
  }

       
  bool Structure :: Load( const TiXmlElement &_element )
  {
    const TiXmlElement *child, *parent;
    double d; atat::rVector3d vec;
    int i;

    // Find first XML "Structure" node (may be _element from start, or a child of _element)
    std::string str = _element.Value();
    if ( str.compare("Structure" ) != 0 )
      parent = _element.FirstChildElement("Structure");
    else
      parent = &_element;
    __DOASSERT( not parent, "Could not find Structure tag in xml.\n" )

    // read PI name if available
    name = "";
    if ( parent->Attribute("name") )
      name = parent->Attribute("name");
    energy = 1.0;
    if ( not parent->Attribute("energy", &d) )
      energy = 666.666;
    energy = types::t_real(d);

    // reads in cell
    child = parent->FirstChildElement( "Cell" );
    __DOASSERT( not child, "Unit-cell not found in structure input\n")
    child = child->FirstChildElement( "row" );
    freeze = FREEZE_NONE;
    for (i=0 ; child and i<3; child=child->NextSiblingElement( "row" ), i++ )
    {
      d=1.0; __DOASSERT( not child->Attribute("x", &d),
                         "Incomplete cell in structure tag.\n" ); vec(0) = d;
      d=1.0; __DOASSERT( not child->Attribute("y", &d),
                         "Incomplete cell in structure tag.\n"); vec(1) = d;
      d=1.0; __DOASSERT( not child->Attribute("z", &d),
                         "Incomplete cell in structure tag.\n"); vec(2) = d;
      cell.set_row(i,vec);
      if ( child->Attribute("freeze") )
      {
        str = child->Attribute("freeze");
        if ( str.find("all") != std::string::npos )
          switch( i )
          {
            default:
            case 0: freeze |= FREEZE_XX | FREEZE_XY | FREEZE_XZ; break;
            case 1: freeze |= FREEZE_XY | FREEZE_YY | FREEZE_YZ; break;
            case 2: freeze |= FREEZE_XZ | FREEZE_YZ | FREEZE_ZZ; break;
          }
        else
        {
          if ( str.find("x") != std::string::npos )
            switch( i )
            {
              default:
              case 0: freeze |= FREEZE_XX; break;
              case 1: freeze |= FREEZE_XY; break;
              case 2: freeze |= FREEZE_XZ; break;
            }
          if ( str.find("y") != std::string::npos )
            switch( i )
            {
              default:
              case 0: freeze |= FREEZE_XY; break;
              case 1: freeze |= FREEZE_YY; break;
              case 2: freeze |= FREEZE_YZ; break;
            }
          if ( str.find("z") != std::string::npos )
            switch( i )
            {
              default:
              case 0: freeze |= FREEZE_XZ; break;
              case 1: freeze |= FREEZE_YZ; break;
              case 2: freeze |= FREEZE_ZZ; break;
            }
        }  
      }
    }
    __DOASSERT(i != 3, "More than three row for unit-cell in input\n")

    scale = 0;
    parent->Attribute("scale", &scale);
    if ( not scale )
    {
      scale = 2.0;
      if ( lattice )
       scale = lattice->scale; 
    }

    // reads in atoms
    child = parent->FirstChildElement( "Atom" );
    atoms.clear();
    for (; child; child=child->NextSiblingElement( "Atom" ) )
    {
      StrAtom stratom;
      Atom atom;
      if( not (     lattice
                and stratom.Load( *child ) 
                and lattice->convert_StrAtom_to_Atom( stratom, atom )  ) )
        __TRYASSERT( not atom.Load( *child ),
                     "Could not read atom from input.\n" )

      if (    ( lattice and atom.site > (types::t_int)lattice->sites.size() )
           or atom.site < -1 )
        atom.site = -1;
      if ( lattice and atom.site < 0 )
        atom.site = lattice->get_atom_site_index( atom );
      if ( lattice and atom.site >= 0 )
        atom.freeze |= lattice->sites[ atom.site ].freeze;
      atoms.push_back(atom);
    }


    // reads in kvectors
    child = parent->FirstChildElement( "Kvec" );
    k_vecs.clear();
    for (; child; child=child->NextSiblingElement( "Kvec" ) )
    {
      CAtom kvec;
      __DOASSERT( not kvec.Load(*child),
                  "Could not load k-vector from xml input.\n" )
      k_vecs.push_back(kvec);
    }
    if ( lattice and ( not k_vecs.size() ) )
      find_k_vectors();

    return true;
  }

  bool Structure :: set_site_indices()
  {
    if ( not lattice )
      return false;

    bool result = true;
    t_Atoms :: iterator i_atom = atoms.begin();
    t_Atoms :: iterator i_atom_end = atoms.end();
    for(; i_atom != i_atom_end; ++i_atom )
    {
      i_atom->site = lattice->get_atom_site_index( i_atom->pos );
      (i_atom->site == -1) ?
        result = false:
        i_atom->freeze |= lattice->sites[ i_atom->site ].freeze;
    }
    return result;
  }

  // prints an xml Structure tag. 
  // It may or may not have been created by a call
  // to Constituent_Stain :: print_xml( ... )
  void Structure :: print_xml( TiXmlElement &_node ) const
  {
    TiXmlElement *child, *parent, *structure;

    // insert cell
    structure = &_node;
    std::string name = structure->Value(); 
    if ( name.compare("Structure") )
    {
      structure = new TiXmlElement( "Structure" );
      _node.LinkEndChild( structure );
    }
    parent = new TiXmlElement( "Cell" );
    structure->LinkEndChild( parent );
    structure->SetAttribute("N", atoms.size() );
    structure->SetAttribute("N", atoms.size() );
    structure->SetAttribute("name", name );
    
    for (int i=0; i < 3; ++i)
    {
      child = new TiXmlElement( "row" );
      child->SetDoubleAttribute( "x", cell.get_row(i)(0) );
      child->SetDoubleAttribute( "y", cell.get_row(i)(1) );
      child->SetDoubleAttribute( "z", cell.get_row(i)(2) );
      parent->LinkEndChild(child);
      if ( i == 0 and
           freeze & FREEZE_XX or 
           freeze & FREEZE_XY or
           freeze & FREEZE_XZ )
      {
         std::ostringstream ss(""); 
         if ( freeze & (FREEZE_XX | FREEZE_XY | FREEZE_XZ ) )
           ss << "all";
         else
         {
           if ( freeze & FREEZE_XX )
             ss << "x";
           if ( freeze & FREEZE_XY )
             ss << "y";
           if ( freeze & FREEZE_XZ )
             ss << "z";
         }
         _node.SetAttribute( "freeze", ss.str().c_str() );
      }
      else if ( i == 1 and
                freeze & FREEZE_XY or 
                freeze & FREEZE_YY or
                freeze & FREEZE_YZ )
      {
         std::ostringstream ss(""); 
         if ( freeze & (FREEZE_XY | FREEZE_YY | FREEZE_YZ ) )
           ss << "all";
         else
         {
           if ( freeze & FREEZE_XY )
             ss << "x";
           if ( freeze & FREEZE_YY )
             ss << "y";
           if ( freeze & FREEZE_YZ )
             ss << "z";
         }
         _node.SetAttribute( "freeze", ss.str().c_str() );
      }
      else if ( i == 2 and
                freeze & FREEZE_XZ or 
                freeze & FREEZE_YZ or
                freeze & FREEZE_ZZ )
      {
         std::ostringstream ss(""); 
         if ( freeze & (FREEZE_XZ | FREEZE_YZ | FREEZE_YZ ) )
           ss << "all";
         else
         {
           if ( freeze & FREEZE_XZ )
             ss << "x";
           if ( freeze & FREEZE_YZ )
             ss << "y";
           if ( freeze & FREEZE_ZZ )
             ss << "z";
         }
         _node.SetAttribute( "freeze", ss.str().c_str() );
      }
    } // end of for (int i=0; i < 3; ++i) over cell vectors

    std::vector<Atom> :: const_iterator i_atom = atoms.begin();
    std::vector<Atom> :: const_iterator i_atom_end = atoms.end();
    for (; i_atom != i_atom_end; ++i_atom) 
    {
      TiXmlElement *xml_atom = new TiXmlElement("Atom");
      StrAtom stratom;
      if ( lattice )
        lattice->convert_Atom_to_StrAtom(*i_atom, stratom);
      stratom.print_xml(xml_atom);
      structure->LinkEndChild(xml_atom);
    }

    if ( k_vecs.size() == 0 )
      return; 

    std::vector<CAtom> :: const_iterator i_kvec = k_vecs.begin();
    std::vector<CAtom> :: const_iterator i_kvec_end = k_vecs.end();
    for (; i_kvec != i_kvec_end; ++i_kvec) 
    {
      TiXmlElement *xml_kvec = new TiXmlElement("Kvec");
      i_kvec->print_xml(xml_kvec);
      structure->LinkEndChild(xml_kvec);
    }

  }

   void Structure :: find_k_vectors()
   {
     if ( not lattice ) return;
   
     atat::rVector3d kvec;
     atat::rMatrix3d k_lat = !( lattice->cell );
     atat::rMatrix3d k_cell = !( cell );
     k_vecs.clear();
  
  
     // A is the basis used to determine "a" first brillouin zone
     atat::rMatrix3d A = (!k_lat) * k_cell;
     
     types::t_int range =  2 * types::t_int( atoms.size() );
     // sets up the n-dimensional iterators
     opt::Ndim_Iterator< types::t_int, std::less_equal<types::t_int> > global_iterator;
     global_iterator.add( -range, range);
     global_iterator.add( -range, range);
     global_iterator.add( -range, range);
     
     // the following loop creates all possible k-vectors,
     // it then refolds them and adds them to the k vector list
     // only one copy of each refolded vector is allowed
     k_vecs.clear();
     do
     {
       // creates vector in A basis
       kvec[0] =  (types::t_real) global_iterator.access(0);
       kvec[1] =  (types::t_real) global_iterator.access(1);
       kvec[2] =  (types::t_real) global_iterator.access(2);
       kvec = A * kvec;
     
       // if any of the coordinates is >= 1, then this is a periodic image
       if (    Fuzzy::geq( kvec(0), 1.0 ) 
            or Fuzzy::geq( kvec(1), 1.0 ) 
            or Fuzzy::geq( kvec(2), 1.0 ) ) continue;
       // if any of the coordinates is < 0, then this is a periodic image
       if (    Fuzzy::le( kvec(0), 0.0 ) 
            or Fuzzy::le( kvec(1), 0.0 ) 
            or Fuzzy::le( kvec(2), 0.0 ) ) continue;
      
       // Goes back to lattice basis
       kvec[0] =  (types::t_real) global_iterator.access(0);
       kvec[1] =  (types::t_real) global_iterator.access(1);
       kvec[2] =  (types::t_real) global_iterator.access(2);
       // And then to cartesian
       kvec = k_cell * kvec;

       // Refold centers all vectors around origin, eg vector within fist
       // brillouin zone
       refold(kvec, k_lat);

       k_vecs.push_back( t_kAtom(kvec,0) );
     
     } while( ++global_iterator );
  
  
     // finally, sorts  k_vec according to size
     std::sort( k_vecs.begin(), k_vecs.end(), sort_kvec );
   }


  void  find_range( const atat::rMatrix3d &A, atat::iVector3d &kvec )
  {
    atat::rVector3d a = A.get_row(0), b;
    types::t_int n = 1;
    b = a;
    while( not is_int(b) )
      { b += a; n++;  }
    kvec[0] = n;

    a = A.get_row(1);
    b = a; n = 1;
    while( not is_int(b) )
      { b += a; n++;  }
    kvec[1] = n;
    
    a = A.get_row(2);
    b = a; n = 1;
    while( not is_int(b) )
      { b += a; n++;  }
    kvec[2] = n;
  }

  bool are_equivalent( const atat::rVector3d &_a,
                       const atat::rVector3d &_b,
                       const atat::rMatrix3d &_cell) 
  {

    opt::Ndim_Iterator<types::t_int, std::less_equal<types::t_int> > i_cell;
    atat::rVector3d compute;

    i_cell.add(-1,1);
    i_cell.add(-1,1);
    i_cell.add(-1,1);

    do
    {
      compute(0) = (types::t_real) i_cell.access(0);
      compute(1) = (types::t_real) i_cell.access(1);
      compute(2) = (types::t_real) i_cell.access(2);

      compute = _a + _cell*compute;
      if ( norm2( compute - _b ) <  types::tolerance ) 
        return true;

    } while ( (++i_cell) );

    return false;
  }

  std::ostream& Structure :: print_xcrysden( std::ostream &_stream ) const
  {
    if( not lattice ) return _stream;
    _stream << "CRYSTAL\nPRIMVEC\n" << ( (~cell) * scale ) << "PRIMCOORD\n" 
              << atoms.size() << " 1 \n";  
    t_Atoms :: const_iterator i_atom = atoms.begin();
    t_Atoms :: const_iterator i_atom_end = atoms.end();
    bool which = true;
    for(; i_atom != i_atom_end; ++i_atom, which = not which )
    {
      StrAtom stratom; lattice->convert_Atom_to_StrAtom( *i_atom, stratom );
      _stream << " " << Physics::Atomic::Z(stratom.type) 
              << " " << ( stratom.pos * scale ) << "\n";
    }
    return _stream;
  } 

  std::ostream& Structure :: print_xyz( std::ostream &_stream,
                                        const std::string &_name ) const
  {
    if( not lattice ) return _stream;
    _stream << atoms.size() << "\n"
            << _name << "\n";
    t_Atoms :: const_iterator i_atom = atoms.begin();
    t_Atoms :: const_iterator i_atom_end = atoms.end();
    bool which = true;
    for(; i_atom != i_atom_end; ++i_atom, which = not which )
    {
      StrAtom stratom; lattice->convert_Atom_to_StrAtom( *i_atom, stratom );
      _stream << "   " << std::setw(2) << stratom.type
              << "   " << ( stratom.pos * scale ) << "\n";
    }
    return _stream;
  } 

} // namespace Crystal


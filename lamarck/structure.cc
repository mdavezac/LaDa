#include <algorithm>
#include <functional>

#ifndef __PGI
  #include<ext/functional>
  using __gnu_cxx::compose1;
#else
  #include<functional>
  using std::compose1;
#endif

#include <opt/compose_functors.h>

#include "atat/misc.h"
#include "structure.h"


namespace Ising_CE {

  using atat::rndseed;
  Lattice* Structure :: lattice = NULL;

  const types::t_unsigned Structure::FREEZE_NONE = 0;
  const types::t_unsigned Structure::FREEZE_XX = 1;
  const types::t_unsigned Structure::FREEZE_XY = 2;
  const types::t_unsigned Structure::FREEZE_XZ = 4;
  const types::t_unsigned Structure::FREEZE_YY = 8;
  const types::t_unsigned Structure::FREEZE_YZ = 16;
  const types::t_unsigned Structure::FREEZE_ZZ = 32;


  void Structure :: convert_from_ATAT ( const Atat_Structure &atat  )
  {
    cell = atat.cell;
    for ( types::t_int i = 0; i < atat.atom_pos.get_size(); i++ )
      atoms.push_back( Atom_Type<types::t_real>( atat.atom_pos[i],
                             static_cast<types::t_real>(atat.atom_type[i]) ) );
  }
                                                                 
  void Structure :: print_out (std::ostream &stream) const
  {
    stream << std::endl << " Structure Cell " << std::endl;
    stream << cell;
    
    #ifdef _DEBUG_LADA_
      if (atoms.size() == 0 )
      {
        std::cerr << " Structure contains no atoms!! " << std::endl;
        exit(0);
      }
    #endif
    stream << " Structure atoms " << std::endl;
    std::vector<Atom> :: const_iterator i_atom = atoms.begin();
    std::vector<Atom> :: const_iterator i_end = atoms.end();
    for( ; i_atom != i_end; ++i_atom )
    {
      stream << "  Position: ";
      i_atom->print_out(stream);
    }
    if ( not k_vecs.size() ) 
      return;
    stream << " Structure K vectors " << std::endl;
    std::vector<CAtom> :: const_iterator i_kvec = k_vecs.begin();
    std::vector<CAtom> :: const_iterator i_kvec_end = k_vecs.end();
    for( ; i_kvec != i_kvec_end; ++i_kvec )
    {
      stream << "  Kvec: ";
      i_kvec->print_out(stream);
    }
  }

  void Structure :: set_atom_types( const std::vector<types::t_real> &types)
  {
    #ifdef _DEBUG_LADA_
      if ( types.size() != atoms.size() )
      {
        std::cerr << "Vectors are not of equivalent size in "
                  << "void set_atom_types( const std::vector<types::t_real> &)"
                  << std::endl;
        exit(0);
      }
    #endif // _DEBUG_LADA_

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
    #ifdef _DEBUG_LADA_
      if ( types.size() != atoms.size() )
      {
        std::cerr << "Vectors are not of equivalent size in "
                  << "void set_atom_types( const std::vector<types::t_real> &)"
                  << std::endl;
        exit(0);
      }
    #endif // _DEBUG_LADA_

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

  void Structure :: randomize (types::t_real range, bool centered = false)
  {
    std::vector<Atom> :: iterator i_atom = atoms.begin();
    std::vector<Atom> :: iterator i_last = atoms.end();
    rndseed();
    if (  centered )
    {
      for(; i_atom != i_last; i_atom++ )
      #ifdef WANG_CONSTRAINTS 
        {
           i_atom->type = 1.0 - (types::t_real) ( rand()/( (types::t_real)RAND_MAX ) ) * range * 0.5;
           if ( rand() / ( (types::t_real)RAND_MAX ) < 0.5 ) 
             i_atom->type *= -1.0;
        }
      #else
        {
           i_atom->type = ( rand() / ( (types::t_real)RAND_MAX ) < 0.5 ) ? -1.0 : 1.0;
           i_atom->type += (types::t_real) ( rand()/( (types::t_real)RAND_MAX ) - 0.5 ) * range;
        }
      #endif // WANG_CONSTRAINTS
    }
    else 
    {
      for(; i_atom != i_last; i_atom++ )
        i_atom->type = (types::t_real) ( rand()/( (types::t_real)RAND_MAX ) -0.5 )
                       * range;
    }
  }
       
  bool Structure :: Load( const TiXmlElement &_element )
  {
    const TiXmlElement *child, *parent;
    double d; atat::rVector3d vec;
    int i;

    std::string str = _element.Value();
    if ( str.compare("Structure" ) != 0 )
      parent = _element.FirstChildElement("Structure");
    else
      parent = &_element;
    if ( not parent )
      return false;

    // read PI name if available
    if ( not parent->Attribute("PI", &i) )
      Pi_name = 0;
    Pi_name = types::t_int(i);
    energy = 1.0;
    if ( not parent->Attribute("energy", &d) )
      energy = 666.666;
    energy = types::t_real(d);

    // reads in cell
    child = parent->FirstChildElement( "Cell" );
    if ( !child )
      return false;
    child = child->FirstChildElement( "column" );
    freeze = FREEZE_NONE;
    for (i=0 ; child and i<3; child=child->NextSiblingElement( "column" ), i++ )
    {
      d=1.0; if( not child->Attribute("x", &d) ) return false; vec(0) = d;
      d=1.0; if( not child->Attribute("y", &d) ) return false; vec(1) = d;
      d=1.0; if( not child->Attribute("z", &d) ) return false; vec(2) = d;
      cell.set_column(i,vec);
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
    if (i != 3)
      return false;

    scale = 0;
    parent->Attribute("scale", &scale);
    if ( not scale )
    {
      scale = 2.0;
      if ( lattice )
       scale = lattice->scale; 
    }

    // reads in atoms
    types::t_real (*ptr_norm)(const atat::FixedVector<types::t_real, 3> &) = &atat::norm2;
    child = parent->FirstChildElement( "Atom" );
    Atom atom;
    StrAtom stratom;
    atoms.clear();
    for (; child; child=child->NextSiblingElement( "Atom" ) )
    {
      if ( not stratom.Load(*child) )
       return false;
      if ( (not lattice) or (not lattice->convert_StrAtom_to_Atom( stratom, atom )) )
        if ( not atom.Load( *child ) )
          return false;
      atoms.push_back(atom);
    }
//   if ( atoms.size() )
//     std::sort( atoms.begin(), atoms.end(), 
//                opt::ref_compose2( std::less<types::t_real>(), std::ptr_fun( ptr_norm ),
//                                   std::ptr_fun( ptr_norm ) ) );

    // reads in atoms
    child = parent->FirstChildElement( "Kvec" );
    k_vecs.clear();
    for (; child; child=child->NextSiblingElement( "Kvec" ) )
    {
      CAtom kvec;
      if ( not kvec.Load(*child) )
        return false;
      k_vecs.push_back(kvec);
    }

    if ( k_vecs.size() )
      std::sort( k_vecs.begin(), k_vecs.end(), 
                 opt::ref_compose2( std::less<types::t_real>(), std::ptr_fun( ptr_norm ),
                                    std::ptr_fun( ptr_norm ) ) );
    return true;
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
    
    for (int i=0; i < 3; ++i)
    {
      child = new TiXmlElement( "column" );
      child->SetDoubleAttribute( "x", cell.get_column(i)(0) );
      child->SetDoubleAttribute( "y", cell.get_column(i)(1) );
      child->SetDoubleAttribute( "z", cell.get_column(i)(2) );
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

// void Structure :: operator=( const Structure &_str )
// {
//   cell = _str.cell;
//   atoms.clear(); atoms.resize( _str.atoms.size() );
//   std::copy( _str.atoms.begin(), _str.atoms.end(), atoms.begin() );
//   k_vecs.clear(); k_vecs.resize( _str.k_vecs.size() );
//   std::copy( _str.k_vecs.begin(), _str.k_vecs.end(), k_vecs.begin() );
//   Pi_name = _str.Pi_name;
//   energy = _str.energy;
// }
} // namespace Ising_CE

#ifdef _MPI
#include<mpi/mpi_object.h>

namespace mpi
{
  template<>
  bool BroadCast :: serialize<Ising_CE::Structure>( Ising_CE::Structure &_str )
  {
    // first copies cell vectors
    if ( not serialize( _str.cell.x[0], (_str.cell.x[0]+3) ) ) return false;
    if ( not serialize( _str.cell.x[1], (_str.cell.x[1]+3) ) ) return false;
    if ( not serialize( _str.cell.x[2], (_str.cell.x[2]+3) ) ) return false;
    
    // copies freeze, Pi_name, energy, and scale
    if ( not serialize( _str.freeze ) ) return false;
    if ( not serialize( _str.Pi_name ) ) return false;
    if ( not serialize( _str.energy ) ) return false;
    if ( not serialize( _str.scale ) ) return false;

    // copies size of atoms and k_vecs
    types::t_int n = _str.atoms.size();
    if( not serialize( n ) ) return false;
    if ( stage == COPYING_FROM_HERE )
      _str.atoms.resize(n);
    n = _str.k_vecs.size();
    if( not serialize( n ) ) return false;
    if ( stage == COPYING_FROM_HERE )
      _str.k_vecs.resize(n);

    // copies atoms
    Ising_CE::Structure::t_Atoms :: iterator i_atom = _str.atoms.begin();
    Ising_CE::Structure::t_Atoms :: iterator i_atom_end = _str.atoms.end();
    for(; i_atom != i_atom_end; ++i_atom )
    {
      if ( not serialize( i_atom->pos.x, i_atom->pos.x+3 ) ) return false;
      if ( not serialize( i_atom->type ) ) return false;
      if ( not serialize( i_atom->freeze ) ) return false;
    }

    // copies k_vecs
    Ising_CE::Structure::t_kAtoms :: iterator i_kvec = _str.k_vecs.begin();
    Ising_CE::Structure::t_kAtoms :: iterator i_kvec_end = _str.k_vecs.end();
    for(; i_kvec != i_kvec_end; ++i_kvec )
    {
      if ( not serialize( i_kvec->pos.x, i_kvec->pos.x+3 ) ) return false;
      if ( not serialize( std::real(i_kvec->type) ) ) return false;
      if ( not serialize( std::imag(i_kvec->type) ) ) return false;
      if ( not serialize( i_kvec->freeze ) ) return false;
    }

    return true;
  }
}

#endif

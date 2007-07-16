#include <algorithm>
#include <functional>

#ifndef __PGI
  #include<ext/functional>
  using __gnu_cxx::compose1;
#else
  #include<functional>
  using std::compose1;
#endif

#include "opt/compose_functors.h"
#include "opt/ndim_iterator.h"

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
      stream << std::endl;
    }
    if ( not k_vecs.size() ) 
      return;
    stream << " Structure K vectors " << std::endl;
    std::vector<CAtom> :: const_iterator i_kvec = k_vecs.begin();
    std::vector<CAtom> :: const_iterator i_kvec_end = k_vecs.end();
    for( ; i_kvec != i_kvec_end; ++i_kvec )
    {
      stream << "  Kvec: " << std::fixed << std::setprecision(5);
      i_kvec->print_out(stream);
      stream << std::endl;
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
      if (    ( lattice and atom.site > (types::t_int)lattice->sites.size() )
           or atom.site < -1 )
        atom.site = -1;
      if ( lattice and atom.site != -1 )
        atom.freeze |= lattice->sites[ atom.site ].freeze;
      atoms.push_back(atom);
    }

    // reads in kvectors
    child = parent->FirstChildElement( "Kvec" );
    k_vecs.clear();
    for (; child; child=child->NextSiblingElement( "Kvec" ) )
    {
      CAtom kvec;
      if ( not kvec.Load(*child) )
        return false;
      k_vecs.push_back(kvec);
    }
    if ( lattice and ( not k_vecs.size() ) )
      find_k_vectors();

    if ( k_vecs.size() )
      std::sort( k_vecs.begin(), k_vecs.end(), 
                 opt::ref_compose2( std::less<types::t_real>(), std::ptr_fun( ptr_norm ),
                                    std::ptr_fun( ptr_norm ) ) );
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

  void Structure :: find_k_vectors()
  {
    if ( not lattice ) return;
  
    atat::rVector3d kvec;
    atat::rMatrix3d k_lat = !( ~(lattice->cell) );
    atat::rMatrix3d k_cell = !( ~(cell) );
    k_vecs.clear();


    // A is the basis used to determine "a" first brillouin zone
    atat::rMatrix3d A = (!k_lat) * k_cell;
    atat::iVector3d range, min;
    
    // computes range up to first periodic image
    find_range( A, range );
    
    // sets up the n-dimensional iterators
    opt::Ndim_Iterator< types::t_int, std::less_equal<types::t_int> > global_iterator;
    global_iterator.add( 0, range[0]);
    global_iterator.add( 0, range[1]);
    global_iterator.add( 0, range[2]);
    
    // the following loop creates all possible k-vectors,
    // it then refolds them and adds them to the k vector list
    // only one copy of each refolded vector is allowed
    types::t_real (*ptr_norm)(const atat::FixedVector<types::t_real, 3> &) = &atat::norm2;
    t_kAtoms :: iterator i_begin = k_vecs.begin();
    t_kAtoms :: iterator i_end = k_vecs.end();
    t_kAtoms :: iterator i_which;
    do
    {
      // creates vector in A basis
      kvec[0] =  (types::t_real) global_iterator.access(0);
      kvec[1] =  (types::t_real) global_iterator.access(1);
      kvec[2] =  (types::t_real) global_iterator.access(2);
      kvec = A * kvec;
    
      kvec[0] -= rint(kvec[0]); 
      kvec[1] -= rint(kvec[1]); 
      kvec[2] -= rint(kvec[2]); 
      
      // switches to cartesian coordinates
      kvec = k_lat * kvec;
      
      // find if vector is already in list
      i_which = std::find_if( i_begin, i_end, 
                     compose1( std::bind2nd(std::less<types::t_real>(), types::tolerance),
                     compose1( std::ptr_fun(ptr_norm),
                               bind2nd(std::minus<atat::rVector3d>(), kvec) ) ) );
      // if it is in list, don't add it
      if ( i_which != i_end  ) 
        continue;
      
      k_vecs.push_back( t_kAtom(kvec,0) );
      i_begin = k_vecs.begin();
      i_end = k_vecs.end();
    
    } while( ++global_iterator );

    // refolds the vectors somewhat better
    i_begin = k_vecs.begin();
    i_end = k_vecs.end();
    for( ; i_begin != i_end; i_begin++ )
      refold(i_begin->pos, k_lat);

    // finally, puts vector 0,0,0 at front of list
    i_begin = k_vecs.begin();
    i_end = k_vecs.end();
    i_which = std::find_if( i_begin, i_end, 
                   compose1( std::bind2nd(std::less<types::t_real>(), types::tolerance),
                   compose1( std::ptr_fun(ptr_norm), std::_Identity<atat::rVector3d>() ) ) );
 
    if ( i_which != i_end  ) 
      std::iter_swap( i_which, k_vecs.begin() );
      
    // the refolding is not perfect, we now remove equivalent
    // vectors "by hand "
    remove_equivalents(k_vecs, k_lat);
  }

  void  find_range( const atat::rMatrix3d &A, atat::iVector3d &kvec )
  {
    atat::rVector3d a = A.get_column(0), b;
    types::t_int n = 1;
    b = a;
    while( not is_int(b) )
      { b += a; n++;  }
    kvec[0] = n;

    a = A.get_column(1);
    b = a; n = 1;
    while( not is_int(b) )
      { b += a; n++;  }
    kvec[1] = n;
    
    a = A.get_column(2);
    b = a; n = 1;
    while( not is_int(b) )
      { b += a; n++;  }
    kvec[2] = n;
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

} // namespace Ising_CE

#ifdef _MPI
#include<mpi/mpi_object.h>

namespace mpi
{
  template<>
  bool BroadCast :: serialize<Ising_CE::Structure>( Ising_CE::Structure &_str )
  {
    // first copies cell vectors
    if ( not serialize( _str.cell ) ) return false;
    
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
      if ( not serialize( i_atom->pos ) ) return false;
      if ( not serialize( i_atom->type ) ) return false;
      if ( not serialize( i_atom->freeze ) ) return false;
    }

    // copies k_vecs
    Ising_CE::Structure::t_kAtoms :: iterator i_kvec = _str.k_vecs.begin();
    Ising_CE::Structure::t_kAtoms :: iterator i_kvec_end = _str.k_vecs.end();
    for(; i_kvec != i_kvec_end; ++i_kvec )
    {
      if ( not serialize( i_kvec->pos ) ) return false;
      if ( not serialize( std::real(i_kvec->type) ) ) return false;
      if ( not serialize( std::imag(i_kvec->type) ) ) return false;
      if ( not serialize( i_kvec->freeze ) ) return false;
    }

    return true;
  }
}

#endif
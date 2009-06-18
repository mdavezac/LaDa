//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <algorithm>
#include <sstream>
#include <boost/lexical_cast.hpp>
#include <boost/tuple/tuple_io.hpp>

#include <opt/ndim_iterator.h>
#include <opt/debug.h>
#include <opt/tinyxml.h>

#include <atat/misc.h>
#include <physics/physics.h>

#include "structure.h"
#include "fill_structure.h"
#include "epi_structure.h"
#include "fourier.h"
#include "smith.h"

#include <opt/smith_normal_form.h>


namespace LaDa
{

  namespace Crystal 
  {

    namespace details
    {
      Crystal::Lattice* Structure :: lattice = NULL;
    }

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
             << "   " << std::setw(9) << cell(0,1)
             << "   " << std::setw(9) << cell(0,2) << "\n"
             << "   " << std::setw(9) << cell(1,0)
             << "   " << std::setw(9) << cell(1,1)
             << "   " << std::setw(9) << cell(1,2) << "\n"
             << "   " << std::setw(9) << cell(2,0)
             << "   " << std::setw(9) << cell(2,1)
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

         
    const TiXmlElement * Structure :: find_node( const TiXmlElement &_element )
    {
      const TiXmlElement *parent = opt::find_node( _element, "Structure" );
      __DOASSERT( not parent, "Could not find Structure tag in xml.\n" )
      return parent;
    }

    void Structure :: load_attributes( const TiXmlElement &_element )
    {
      // read PI name if available
      name = "";
      if ( _element.Attribute("name") )
        name = _element.Attribute("name");
      energy = 666.666e0;
      if ( _element.Attribute("energy") )
        energy = boost::lexical_cast<types::t_real>( _element.Attribute("energy") );
      if ( _element.Attribute("weight") )
        weight = boost::lexical_cast<types::t_real>( _element.Attribute("weight") );

      scale = 2.0;
      if ( _element.Attribute("scale") )
        scale = boost::lexical_cast<types::t_real>( _element.Attribute("scale") );
      else if ( lattice )
         scale = lattice->scale; 
    }

    bool Structure :: load_cell( const TiXmlElement &_element )
    {
      // reads in cell
      const TiXmlElement* child = _element.FirstChildElement( "Cell" );
      if( not child ) return false;
      child = child->FirstChildElement( "row" );
      freeze = FREEZE_NONE;
      size_t i(0);
      for (; child and i<3; child=child->NextSiblingElement( "row" ), i++ )
      {
        __DOASSERT( not (     child->Attribute("x") 
                          and child->Attribute("y") 
                          and child->Attribute("z")  ),
                    "Incomplete cell in xml input.\n" )

        __TRYBEGIN
          cell(i,0) = boost::lexical_cast<types::t_real>( child->Attribute("x") );
          cell(i,1) = boost::lexical_cast<types::t_real>( child->Attribute("y") );
          cell(i,2) = boost::lexical_cast<types::t_real>( child->Attribute("z") );
        __TRYEND(, "Could not parse cell input.\n" )
        if ( child->Attribute("freeze") )
        {
          std::string str = child->Attribute("freeze");
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
      return i == 3;
    }

    bool Structure :: load_atoms( const TiXmlElement &_element )
    {
      // reads in atoms
      const TiXmlElement* child = _element.FirstChildElement( "Atom" );
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
      return atoms.size() != 0;
    }

    bool Structure :: load_kvecs( const TiXmlElement &_element )
    {
      // reads in kvectors
      const TiXmlElement* child = _element.FirstChildElement( "Kvec" );
      k_vecs.clear();
      for (; child; child=child->NextSiblingElement( "Kvec" ) )
      {
        CAtom kvec;
        __DOASSERT( not kvec.Load(*child),
                    "Could not load k-vector from xml input.\n" )
        k_vecs.push_back(kvec);
      }

      return k_vecs.size() != 0;
    }

    bool Structure :: load_epitaxial( const TiXmlElement &_node )
    {
      namespace bt = boost::tuples;
      if(     ( not _node.Attribute("direction") )
          and ( not _node.Attribute("extent") ) ) return false;
      if(     ( not _node.Attribute("direction") )
          or  ( not _node.Attribute("extent") ) )
      {
        std::cerr << "Either cell direction, multiplicity, "
                     "or layersize are missing on input\n";
        return false;
      }

      atat::rVector3d direction;
      atat::iVector3d extent;
      
      // First, Load Attributes 
      __TRYBEGIN
        std::istringstream sstr; sstr.str( _node.Attribute("direction") );
        boost::tuple< types::t_real&, types::t_real&, types::t_real& >
                   tupledir( direction[0], direction[1], direction[2] );
        sstr >> boost::tuples::set_open('(')
             >> boost::tuples::set_close(')')
             >> tupledir;
        __DOASSERT( atat::norm2( direction ) < types::tolerance,
                    "direction cannot be null.\n" )
        sstr.str( _node.Attribute("extent") );
        boost::tuple< types::t_int&, types::t_int&, types::t_int& >
                   tupleext( extent[0], extent[1], extent[2] );
        sstr >> boost::tuples::set_open('(')
             >> boost::tuples::set_close(')')
             >> tupleext;
        __DOASSERT( atat::norm2( direction ) < types::tolerance,
                    "extent cannot be null.\n" )
    
        direction = (!lattice->cell) * direction;
        
      __TRYEND(,"Error while loading superlattice description.\n" )
    
      return create_epitaxial_structure( *this, direction, extent );
    }

    bool Structure :: Load( const TiXmlElement &_element )
    {
      const TiXmlElement* parent = find_node( _element );
      __DOASSERT( not parent, "Could not find structure tag.\n" )

      if( parent->Attribute("filename") )
      {
        namespace bfs = boost::filesystem;
        const bfs::path path( Print::reformat_home( parent->Attribute( "filename" ) ) );
        __DOASSERT( not bfs::exists( path ), path.string() + " does not exist.\n" )
        TiXmlDocument doc;
        opt::read_xmlfile( path, doc );
        __DOASSERT( not doc.FirstChild( "Job" ),
                    "Root tag <Job> does not exist in " + path.string() + ".\n" )
        parent = opt::find_node( *doc.FirstChildElement( "Job" ),
                                 "Structure" );
      
        if( parent ) return Load( *parent );
        std::cerr << "Could not find an <Structure> tag in input file.\n";
        return false;
      }
      
      
      // trie to load cell.
      if( load_cell( *parent ) )
      {
        // tries to load atoms.
        bool atomsloaded =  load_atoms( *parent );
        if( not atomsloaded ) atomsloaded = fill_structure( *this );
        __DOASSERT(  not atomsloaded, "Could not read or create atoms.\n" )
      }
      else if( not load_epitaxial( *parent ) ) return false;

      load_attributes( *parent );

      if( not load_kvecs( *parent ) ) find_k_vectors();

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
      std::string nodename = structure->Value(); 
      if ( nodename.compare("Structure") )
      {
        structure = new TiXmlElement( "Structure" );
        _node.LinkEndChild( structure );
      }
      parent = new TiXmlElement( "Cell" );
      structure->LinkEndChild( parent );
      structure->SetAttribute("N", atoms.size() );
      if( Fuzzy::neq( energy, 666.666 ) ) 
        structure->SetDoubleAttribute("energy", energy );
      if( Fuzzy::neq( weight, 1e0 ) )
        structure->SetDoubleAttribute("weight", weight );
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
       namespace bt = boost::tuples;
       if ( not lattice ) return;
      
       k_vecs.clear();
       atat::rMatrix3d const kcell( !(~cell) );
       atat::rMatrix3d const klat( !(~lattice->cell) );
       t_SmithTransform transform = get_smith_transform( kcell, klat );
       
       atat::iVector3d &smith = bt::get<1>(transform);
       const atat::rMatrix3d factor
       ( 
          (~lattice->cell) * (!bt::get<0>(transform))
       );
       for( size_t i(0); i < smith(0); ++i )
         for( size_t j(0); j < smith(1); ++j )
           for( size_t k(0); k < smith(2); ++k )
           {
             // in supercell fractional
             const atat::rVector3d vec1( factor * atat::rVector3d(i,j,k) );
             // in supercell fractional and c entered.
             const atat::rVector3d vec2       
             (                                
               vec1(0) - std::floor( vec1(0) + 0.5 ),
               vec1(1) - std::floor( vec1(1) + 0.5 ),
               vec1(2) - std::floor( vec1(2) + 0.5 )
             );
             // in cartesian
             atat::rVector3d vec( klat * vec2 );
             refold(vec, klat);
           
             k_vecs.push_back( t_kAtom(vec,0) );
           }
     
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
      _stream << "CRYSTAL\nPRIMVEC\n" << ( (~cell) * scale ) << "\nPRIMCOORD\n" 
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

    bool Structure :: set_site_indices()
    {
      if ( not lattice ) return false;

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

    types::t_real concentration( const Structure& _structure, const size_t i )
    {
      __ASSERT( not _structure.lattice, "Lattice pointer is not set.\n" )
      __ASSERT( i >= _structure.lattice->sites.size(), "Index out-of-range.\n" )
      __ASSERT(     _structure.lattice->sites[i].type.size() != 2
                and _structure.lattice->sites[i].type.size() != 1,
                "Wrong number of types at lattice-site.\n" )

      if( _structure.lattice->sites[i].type.size() == 1 ) return 1e0;
      types::t_real result( 0 );
      foreach( const Structure::t_Atom& _a, _structure.atoms )
      {
        __ASSERT( _a.site != 0 and _a.site != 1,
                  "site index not set in structure.\n" )
        if( _a.site == i ) result += _a.type;
      }
      return result;
    }

    void convert_real_to_string_structure( Structure const& _real,
                                           TStructure<std::string> &_string )
    {
      _string.scale = _real.scale;
      _string.freeze = _real.freeze;
      _string.lattice = _real.lattice;
      _string.name = _real.name;
      _string.energy = _real.energy;
      _string.weight = _real.weight;
      _string.cell = _real.cell;
      _string.atoms.resize( _real.atoms.size() );
      TStructure<std::string>::t_Atoms:: iterator i_string = _string.atoms.begin();
      TStructure<std::string>::t_Atoms:: iterator const i_string_end = _string.atoms.end();
      Structure::t_Atoms:: const_iterator i_real = _real.atoms.begin();
      for(; i_string != i_string_end; ++i_string, ++i_real )
        _real.lattice->convert_Atom_to_StrAtom( *i_real, *i_string );
    }

    void convert_string_to_real_structure( TStructure<std::string> const &_string,
                                           Structure &_real )
    {
      _real.scale   = _string.scale;
      _real.freeze  = _string.freeze;
      _real.lattice = _string.lattice;
      _real.name    = _string.name;
      _real.energy  = _string.energy;
      _real.weight  = _string.weight;
      _real.cell    = _string.cell;
      _real.atoms.resize( _string.atoms.size() );
      Structure::t_Atoms:: iterator i_real = _real.atoms.begin();
      Structure::t_Atoms:: iterator const i_real_end = _real.atoms.end();
      TStructure<std::string>::t_Atoms:: const_iterator i_string = _string.atoms.begin();
      for(; i_real != i_real_end; ++i_string, ++i_real )
        _string.lattice->convert_StrAtom_to_Atom( *i_string, *i_real );
      _real.find_k_vectors();
      Fourier::Fourier( _real.atoms.begin(), _real.atoms.end(),
                        _real.k_vecs.begin(), _real.k_vecs.end() );
    }

  } // namespace Crystal

} // namespace LaDa

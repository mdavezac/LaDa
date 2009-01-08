//
//  Version: $Id$
//

#include <boost/lambda/lambda.hpp>
#include <boost/serialization/base_object.hpp>

namespace LaDa
{
  namespace Crystal
  {

    template<class T_TYPE> template< class ARCHIVE >
      void TStructure<T_TYPE> :: serialize( ARCHIVE & _ar, const unsigned int _version)
      {
        _ar & cell;
        _ar & atoms;
        _ar & energy;
        _ar & weight;
        _ar & freeze;
        _ar & scale;
      }
    template< class ARCHIVE >
      void Structure :: serialize( ARCHIVE & _ar, const unsigned int _version)
      {
        _ar & boost::serialization::base_object< TStructure<types::t_real> >(*this);
        _ar & k_vecs;
      }

    //! \cond
    void  find_range( const atat::rMatrix3d &A, atat::iVector3d &kvec );

    template <class CONTAINER>
    void remove_equivalents( CONTAINER &_cont, const atat::rMatrix3d &_cell)
    {
      typename CONTAINER :: iterator i_vec = _cont.begin();
      typename CONTAINER :: iterator i_end = _cont.end();
      typename CONTAINER :: iterator i_which;

      while( i_vec != i_end )
      {
        i_which = i_vec+1;
        for ( ; i_which != i_end; i_which++ )
          if ( are_equivalent( *i_which, *i_vec, _cell ) )
            break;

        if ( i_which == i_end )
        { 
          ++i_vec;
          continue;
        }
        
        ( atat::norm2( (atat::rVector3d&) *i_vec ) < atat::norm2( (atat::rVector3d&) *i_which ) ) ? 
              _cont.erase(i_which): _cont.erase(i_vec);
        i_vec = _cont.begin();
        i_end = _cont.end();
      }

    }
    //! \endcond


    template< class T_TYPE >
     bool TStructure<T_TYPE> :: operator== (const TStructure<T_TYPE> &_str ) const
     {
       return     Fuzzy :: eq( cell(0,0), _str.cell(0,0) )
              and Fuzzy :: eq( cell(1,0), _str.cell(1,0) )
              and Fuzzy :: eq( cell(2,0), _str.cell(2,0) )
              and Fuzzy :: eq( cell(0,1), _str.cell(0,1) )
              and Fuzzy :: eq( cell(1,1), _str.cell(1,1) )
              and Fuzzy :: eq( cell(2,1), _str.cell(2,1) )
              and Fuzzy :: eq( cell(0,2), _str.cell(0,2) )
              and Fuzzy :: eq( cell(1,2), _str.cell(1,2) )
              and Fuzzy :: eq( cell(2,2), _str.cell(2,2) )
              and atoms == _str.atoms;
     }
      
    template<class t_container >
      void Structure :: set_kvectors( const t_container &_container )
      {
        typename t_container :: const_iterator i_kvec =  _container.begin();
        typename t_container :: const_iterator i_end =  _container.end();
        CAtom kvec;
        k_vecs.clear();
        k_vecs.reserve( _container.size() );
        for( ; i_kvec != i_end; ++i_kvec ) 
        {
          kvec = (*i_kvec);
          k_vecs.push_back( kvec );
        }
      }


    template<class T_R_IT, class T_K_IT>
    Fourier :: Fourier( T_R_IT _rfirst, T_R_IT _rend,
                        T_K_IT _kfirst, T_K_IT _kend )
    {
      const std::complex<types::t_real>
         imath(0, -2*3.1415926535897932384626433832795028841971693993751058208);
      
      for (; _kfirst != _kend; ++_kfirst)
      {
        _kfirst->type = std::complex<types::t_real>(0);
        for(T_R_IT i_r( _rfirst ); i_r != _rend; ++i_r )
        {
          _kfirst->type +=   exp( imath * ( i_r->pos[0] * _kfirst->pos[0] +
                                            i_r->pos[1] * _kfirst->pos[1] +
                                            i_r->pos[2] * _kfirst->pos[2] ) )
                           * i_r->type;
        }
      }
    }
    template<class T_R_IT, class T_K_IT, class T_O_IT >
    Fourier :: Fourier( T_R_IT _rfirst, T_R_IT _rend,
                        T_K_IT _kfirst, T_K_IT _kend,
                        T_O_IT _rout ) // sets rvector values from kspace values
    {
      const types::t_complex
         imath(0, 2*3.1415926535897932384626433832795028841971693993751058208);
      for (; _rfirst != _rend; ++_rfirst, ++_rout)
      {
        *_rout = types::t_complex(0,0);
        for(T_K_IT i_k=_kfirst; i_k != _kend; ++i_k)
        {
          *_rout +=   exp( imath * ( _rfirst->pos[0] * i_k->pos[0] +
                                     _rfirst->pos[1] * i_k->pos[1] +
                                     _rfirst->pos[2] * i_k->pos[2] ) )
                    * i_k->type;
        }
      }
    }

    template< class T_TYPE >
      void TStructure<T_TYPE> :: print_out (std::ostream &stream) const
      {
        namespace bl = boost::lambda;
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
                
        
        stream << " Structure atoms\n";
        std::for_each( atoms.begin(), atoms.end(),
                       bl::var(stream) << bl::_1 << "\n" );
      }

    template< class T_TYPE >
      bool TStructure<T_TYPE> :: Load( const TiXmlElement &_element )
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
        weight = 1.0;
        if ( parent->Attribute("energy", &d) )
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
        if ( not scale ) scale = 2.0;
   
        // reads in atoms
        child = parent->FirstChildElement( "Atom" );
        atoms.clear();
        for (; child; child=child->NextSiblingElement( "Atom" ) )
        {
          t_Atom atom;
          __TRYASSERT( atom.Load( *child ), "Could not read atom from input.\n" )
          atoms.push_back(atom);
        }
   
        return true;
      }
      
    template< class T_TYPE >
      void TStructure< T_TYPE > :: print_xml( TiXmlElement &_node ) const
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
        structure->SetDoubleAttribute("energy", energy );
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
   
        typename t_Atoms :: const_iterator i_atom = atoms.begin();
        typename t_Atoms :: const_iterator i_atom_end = atoms.end();
        for (; i_atom != i_atom_end; ++i_atom) 
        {
          TiXmlElement *xml_atom = new TiXmlElement("Atom");
          i_atom->print_xml(xml_atom);
          structure->LinkEndChild(xml_atom);
        }
   
      }

    template< class T_FUNCTIONAL >
      void enumerate_pifile( const std::string &_file, T_FUNCTIONAL &_op )
      {
        __DEBUGTRYBEGIN
        Crystal :: Structure structure;
        std::ifstream file( _file.c_str(), std::ifstream::in );
        do
        {
          if( not Crystal :: read_pifile_structure( file, structure ) ) continue;
          std::cout << "    @" << structure.name << " " 
                    << structure.get_concentration()
                    << " " << _op( structure ) << "\n";
          foreach( Crystal::Structure::t_Atom &atom, structure.atoms )
            atom.type = Fuzzy::gt( atom.type, 0e0 ) ? -1e0: 1e0;
          std::cout << "   -@" << structure.name << " " 
                    << structure.get_concentration()
                    << " " << _op( structure ) << "\n";
        }
        while( not file.eof() );
        __DEBUGTRYEND(, "Error while enumerating pifile.\n" )
      }

  }
} // namespace LaDa

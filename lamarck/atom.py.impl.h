//
//  Version: $Id$
//

  template< class T_TYPE >
  std::string print( const T_TYPE &_at )
  { 
    std::ostringstream sstr;
    _at.print_out( sstr );
    return sstr.str();
  }

  t_Atom* AtomFromObject( list& _ob )
  {
    t_Atom *result = NULL;
    try
    { 
      result = new t_Atom;
      types::t_unsigned length = len(_ob);
      if( length < 3 ) return result;

      result->pos.x[0] = extract< types::t_real >( _ob[0] );
      result->pos.x[1] = extract< types::t_real >( _ob[1] );
      result->pos.x[2] = extract< types::t_real >( _ob[2] );
      if( length == 3 ) return result;
      try{ result->type = extract< types::t_real >(_ob[3] ); }
      catch(...)
      {
        if( not t_Structure::lattice )
          throw std::runtime_error( "Did you forget to initialize the Lattice?" );
        Ising_CE::StrAtom stratom;
        stratom.pos = result->pos;
        stratom.type = extract< std::string >( _ob[3] );
        t_Structure::lattice->convert_StrAtom_to_Atom( stratom, *result );
      }
      if( length == 4 ) return result;
      result->site = extract< types::t_real >( _ob[4] );
      return result;
    }
    catch( std::exception &_e )
    {
      if( result ) delete result;
      std::ostringstream sstr;
      sstr << "Object cannot be converted to an atom: \n"
           << _e.what() << "\n";
      throw std::runtime_error( sstr.str() );
    }
    catch( ... )
    {
      if( result ) delete result;
      throw std::runtime_error( "Could not convert object to Atom." );
    }
    return NULL;
  }
  t_Site* SiteFromObject( list& _ob )
  {
    types::t_unsigned length = len( _ob );
    if( length < 3 )
      throw std::runtime_error( "Object cannot be converted to an atom" );

    t_Site *result = new t_Site;
    try
    { 
      result->pos.x[0] = extract< types::t_real >( _ob[0] );
      result->pos.x[1] = extract< types::t_real >( _ob[1] );
      result->pos.x[2] = extract< types::t_real >( _ob[2] );
      if( length == 3 ) return result;
      list strlist(_ob[3]);
      while( len( strlist ) )
        result->type.push_back( extract<std::string>(strlist.pop(0)) );
      if( length == 4 ) return result;
      result->site = extract< types::t_real >( _ob[4] );
      return result;
    }
    catch( std::exception &_e )
    {
      delete result;
      std::ostringstream sstr;
      sstr << "Object cannot be converted to an atom: \n"
           << _e.what() << "\n";
      throw std::runtime_error( sstr.str() );
    }
    catch( ... )
    {
      delete result;
      throw std::runtime_error( "Could not convert object to Atom." );
    }
    return NULL;
  }



  namespace XML
  {
    typedef Builder< ConstituentStrain::Harmonic::Cubic > t_CubicBuilder;
    typedef t_CubicBuilder :: t_CS t_CubicCS;
    typedef Builder< ConstituentStrain::Harmonic::Tetragonal > t_TetraBuilder;
    typedef t_TetraBuilder :: t_CS t_TetraCS;

    template<> std::string nodename<t_CubicCS>()    { return "CS"; }
    template<> std::string nodename<t_TetraCS>()    { return "CS"; }
    template<> std::string nodename<t_Lattice>()    { return "Lattice"; }
    template<> std::string nodename<t_Structure>()  { return "Structure"; }
    
    template<> void do_specialcode< t_CubicBuilder >( t_CubicBuilder &_type )
      { _type.add_equivalent_clusters(); }
    template<> void do_specialcode< t_TetraBuilder >( t_TetraBuilder &_type )
      { _type.add_equivalent_clusters(); }
    template<> void do_specialcode< t_Lattice >( t_Lattice &_type )
      { t_Structure::lattice = &_type; }
    
    template<> bool doloadcode<t_CubicCS>( t_CubicCS &_type, TiXmlElement *_parent )
      { return _type.Load_Harmonics( *_parent ); }
    template<> bool doloadcode<t_TetraCS>( t_TetraCS &_type, TiXmlElement *_parent )
      { return _type.Load_Harmonics( *_parent ); }
      
    template<> TiXmlElement *findnode<t_CubicCS>( TiXmlHandle &_doc )
    {
      TiXmlElement *result = _doc.FirstChild("Job")
                                 .FirstChild(nodename<t_CubicCS>()).Element();
      if( result ) return result;
      
      result = _doc.FirstChild("Job")
                   .FirstChild("Functional").Element();
      for(; result; result = result->NextSiblingElement("Functional") )
      {
        if( not result->Attribute("type") ) continue;
      
        std::string strg = result->Attribute("type");
        if( strg.compare("CE") == 0 ) break;
      }   
      if( not result ) return NULL;
      result = result->FirstChildElement(nodename<t_CubicCS>());
      return result;
    }
    template<> TiXmlElement *findnode<t_TetraCS>( TiXmlHandle &_doc )
      { return findnode<t_CubicCS>( _doc ); }
    template<> TiXmlElement *findnode<t_CubicBuilder>( TiXmlHandle &_doc )
      { return _doc.FirstChild("Job").Element(); }
    template<> TiXmlElement *findnode<t_TetraBuilder>( TiXmlHandle &_doc )
      { return findnode<t_CubicBuilder>( _doc ); }
  }


  types::t_real toReal(std::string _str )
  { 
    if( not Structure::lattice )
      throw std::runtime_error( "Could not convert atom type.\n" ); 
    if( _str.compare( Structure::lattice->sites[0].type[0] ) == 0 )
      return types::t_real(-1);
    else if( _str.compare( Structure::lattice->sites[0].type[1] ) == 0 )
      return types::t_real(1);
    else
      throw std::runtime_error( "Requested Atomic type is not within lattice" );
  }
  std::string toType( types::t_real _r )
  { 
    if( not Structure::lattice )
      throw std::runtime_error( "Could not convert atom type.\n" ); 
    if( Structure::lattice->sites.size() != 1)
      throw std::runtime_error( "Lattice cannot have more than one atomic site.\n" ); 
    if( Structure::lattice->sites[0].type.size() != 2)
      throw std::runtime_error( "Lattice must have two atomic types.\n" ); 
    if(     Fuzzy::neq( _r, types::t_real(1) ) 
        and Fuzzy::neq( _r, types::t_real(-1) ) ) return "";
    return Fuzzy::neq( _r, types::t_real(1) ) ? 
              Structure::lattice->sites[0].type[0]:
              Structure::lattice->sites[0].type[1] ;
  }

  template< class T_HARMONIC > 
  CEFunc<T_HARMONIC>* generateCEs( Builder<T_HARMONIC> &_builder, t_Structure &_str )
  {
    typedef Builder<T_HARMONIC> t_Builder;
    typedef std::pair< typename t_Builder::t_Chemical*, 
                       typename t_Builder::t_CS*> t_Pair;
    t_Pair pair( _builder.generate_functional(_str) );
    if( (not pair.first) or (not pair.second) ) 
    {
      if( pair.first ) delete pair.first;
      if( pair.second ) delete pair.second;
      throw std::runtime_error( "Could not create functional" );
    }
    CEFunc<T_HARMONIC> *result = new CEFunc<T_HARMONIC>( pair.first, pair.second );
    result->resize( _str.atoms.size() );
    
    return result;
  }

  template< class T_FUNC >
  void assign( T_FUNC &_func, t_Structure &_str )
  {
    std::transform(
       _str.atoms.begin(), _str.atoms.end(), _func.begin(),
       boost::lambda::bind( &t_Structure::t_Atom::type, boost::lambda::_1 ) 
    );
  }

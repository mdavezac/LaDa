//
//  Version: $Id$
//
#ifdef _TYPE_
#  undef _TYPE_
#endif
#ifdef _DIM_
#  undef _DIM_
#endif
#ifdef _CLASSNAME_
#  undef _SHORT_TYPE_ 
#endif
#ifdef _PYTHONNAME_
#  undef _PYTHONNAME_ 
#endif
#ifdef _CONCAT_
#  undef _CONCAT_
#endif
#ifdef _CLASSNAME_
#  undef _CLASSNAME_
#endif

#define _CONCAT_( a, b, c, d ) a ## b ## c ## d

#ifndef _WASINCLUDED_
#  define _WASINCLUDED_ 0
#else
#  undef _WASINCLUDED_
#endif

#if defined(_WASINCLUDED_)

#ifndef _INMODULE_

namespace Ising_CE { namespace details { class PiStructure; } }

typedef Ising_CE::Lattice     t_Lattice;
typedef Ising_CE::Structure   t_Structure;
typedef t_Structure::t_Atom   t_Atom;
typedef t_Lattice::t_Site     t_Site;
#if defined( _CUBIC_CE_ )
  typedef Ising_CE::ConstituentStrain::Harmonic::Cubic      t_Harmonic;
#elif defined( _TETRAGONAL_CE_ )
  typedef Ising_CE::ConstituentStrain::Harmonic::Tetragonal t_Harmonic;
#else 
# error Please define either _CUBIC_CE_ or _TETRAGONAL_CE_
#endif
typedef Ising_CE::ConstituentStrain::Functional<t_Harmonic> t_CS;
typedef VA_CE::Builder<t_Harmonic> t_Builder;
typedef Ising_CE::details::PiStructure t_PiStructure;
typedef opt::ConvexHull::Base<t_PiStructure> t_CH;


namespace Ising_CE
{
  namespace details
  {
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

    template< class T_TYPE > std::string nodename() { return "Unknown Type"; }
    template<> std::string nodename<t_Lattice>()    { return "Lattice"; }
    template<> std::string nodename<t_Structure>()  { return "Structure"; }
    template< class T_TYPE > void do_specialcode( T_TYPE &_type ) {}
    template<> void do_specialcode< t_Builder >( t_Builder &_type )
      { _type.add_equivalent_clusters(); }
    template<> void do_specialcode< t_Lattice >( t_Lattice &_type )
      { t_Structure::lattice = &_type; }
    template< class T_TYPE > bool doloadcode( T_TYPE &_type, TiXmlElement *_parent )
      { return _type.Load( *_parent ); }
    template<> bool doloadcode<t_CS>( t_CS &_type, TiXmlElement *_parent )
      { return _type.Load_Harmonics( *_parent ); }

    template< class T_TYPE > TiXmlElement* findnode( TiXmlHandle &_doc )
      { return _doc.FirstChild("Job").FirstChild(nodename<T_TYPE>()).Element(); }
    template<> TiXmlElement *findnode<t_CS>( TiXmlHandle &_doc )
    {
      TiXmlElement *result = _doc.FirstChild("Job")
                                 .FirstChild(nodename<t_CS>()).Element();
      if( result ) return result;
      
      result = _doc.FirstChild("Job")
                   .FirstChild("functional").Element();
      for(; result; result = result->NextSiblingElement("functional") )
      {
        if( not result->Attribute("type") ) continue;
      
        std::string strg = result->Attribute("type");
        if( strg.compare("CE") == 0 ) break;
      }   
      if( not result ) return NULL;
      result = result->FirstChildElement(nodename<t_CS>());
      return result;
    }
    template<> TiXmlElement *findnode<t_Builder>( TiXmlHandle &_doc )
      { return _doc.FirstChild("Job").Element(); }


    template< class T_TYPE >
    void fromXML(T_TYPE &_type, const std::string &_filename )
    {
      TiXmlDocument doc( _filename ); 
      TiXmlHandle docHandle( &doc ); 

      __DOASSERT( not doc.LoadFile(), 
                     doc.ErrorDesc() << "\n"  
                  << "Could not load input file " << _filename  
                  << ".\nAborting.\n" ) 
      __DOASSERT( not docHandle.FirstChild("Job").Element(),
                  "Could not find <Job> tag in " << _filename << ".\n" )
      TiXmlElement *parent = findnode<T_TYPE>( docHandle );
      __DOASSERT( not parent,    "Could not find <" << nodename<T_TYPE>() 
                              << "> tag in " << _filename << ".\n"   )

      __DOASSERT( not doloadcode( _type, parent ), 
                     "Could not load " << nodename<T_TYPE>()
                  << " from " << _filename << ".\n" )
      do_specialcode( _type );
    }


    void CreateCSfromStructure( t_CS &_cs, t_Structure &_str )
    {
      __DOASSERT( _str.atoms.size() < 1, 
                  "Cannot create constituent strain from "
                  "structure with no atoms.\n" )
      _str.find_k_vectors();
      _cs << _str;
    }

    struct PiStructure
    {
       types::t_int index;
       types::t_real x;
       PiStructure  ( types::t_int _pi = 0, types::t_real _x = -2.0 )
                   : index( _pi ), x( _x ) {}
       types::t_real get_concentration() const { return x; }
       std::string print() const
       { 
         std::ostringstream sstr; 
         sstr << index;
         return sstr.str();
       }
    };
    std::ostream& operator<< ( std::ostream &_os, const PiStructure &_ch )
      { return _os << _ch.index; }

    template< class T_TYPE >
    types::t_unsigned length( const T_TYPE& _type )
      { return _type.size(); }

    class CEFunc : public t_Builder :: t_VA_Functional 
    {
      typedef t_Builder :: t_Chemical t_Chemical;
      typedef t_Builder :: t_CS t_CS;
      typedef t_Builder :: t_VA_functional t_Base;

      public:
        CEFunc() : t_Base() {} 
        CEFunc( t_Chemical* _chem, t_CS* _cs ) : t_Base( _chem, _cs ) {}
        ~CEfunc ()
        { 
          if( _chem ) delete _chem;
          if( _cs ) delete _cs; 
        } 
        types::t_real getitem( types::t_int _i, object ) const;
        {
          if( _i < types::t_int(0) ) _i = size() - _i;
          if( _i >= (types::t_int) size() ) throw std::out_of_range( "CE functional." );
          return (*variables)[ _i ];
        }
        void setitem( types::t_int _i, const object &_object )
        {
          if( _i < types::t_int(0) ) _i = size() - _i;
          if( _i >= (types::t_int) size() ) throw std::out_of_range( "CE functional." );
          types::t_real &var = (*variables)[_i];
          try { var = extract< typename T_TYPE :: t_Type>( _object ); }
          catch(...) 
          {
            if( not Structure::lattice )
              throw std::runtime_error( "Could not convert atom type.\n" ); 
            std::string string =  extract< std::string >( _object );
            if( string.compare( Structure::lattice->sites[0].type[0] ) == 0 )
              var = typename T_TYPE::t_Type(-1);
            else if( string.compare( Structure::lattice->sites[0].type[1] ) == 0 )
              var = typename T_TYPE::t_Type(1);
            else
              throw std::runtime_error( "Requested Atomic type is not within lattice" );
          }
        }
    };

    CEFunc* generateCEs( t_Builder &_builder, t_Structure &_str )
    {
      typedef std::pair<t_Builder::t_Chemical*, t_Builder::t_CS*> t_Pair;
      t_Pair pair( _builder.generate_functional(_str) );
      if( (not pair.first) or (not pair.second) ) 
      {
        if( pair.first ) delete pair.first;
        if( pair.second ) delete pair.second;
        throw std::runtime_error( "Could not create functional" );
      }
      CEFunc *result = new CEFunc( pair.first, pair.second );
      result->resize( _str.atoms.size() );
      
      return result;
    }
  }

}


#else

   class_< t_Structure::t_Atoms >("VecStrings")
     .def(vector_indexing_suite< t_Site::t_Type >());

   class_< t_Atom >( "details_Atom" )
     .def( init< t_Atom >() )
     .def_readwrite( "pos",    &t_Atom::pos )
     .def_readwrite( "site",   &t_Atom::site )
     .def_readwrite( "type",   &t_Atom::type )
     .def_readwrite( "freeze", &t_Atom::freeze )
     .def( "__str__",  &Ising_CE::details::print<t_Atom> ) ;
   class_< t_Site >( "details_Site" )
     .def( init< t_Site >() )
     .def_readwrite( "pos",    &t_Site::pos )
     .def_readwrite( "site",   &t_Site::site )
     .def_readwrite( "type",   &t_Site::type )
     .def_readwrite( "freeze", &t_Site::freeze )
     .def( "__str__",  &Ising_CE::details::print<t_Site> ) ;


   def( "Atom", &Ising_CE::details::AtomFromObject,
        return_value_policy<manage_new_object>() );
   def( "Site", &Ising_CE::details::SiteFromObject,
        return_value_policy<manage_new_object>() );

   class_< t_Structure::t_Atoms >("Atoms")
     .def(vector_indexing_suite< t_Structure::t_Atoms >());
   class_< t_Structure::t_Atoms >("Sites")
     .def(vector_indexing_suite< t_Lattice::t_Sites >());

   class_< t_Structure >( "Structure" )
     .def( init< t_Structure >() )
     .def_readwrite( "cell",   &t_Structure::cell )
     .def_readwrite( "atoms",  &t_Structure::atoms )
     .def_readwrite( "energy", &t_Structure::energy )
     .def_readwrite( "scale",  &t_Structure::scale )
     .def( "__str__",  &Ising_CE::details::print<t_Structure> ) 
     .def( "fromXML",  &Ising_CE::details::fromXML<t_Structure> );

   class_< t_Lattice >( "Lattice" )
     .def( init< t_Lattice >() )
     .def_readwrite( "cell",  &t_Lattice::cell )
     .def_readwrite( "sites", &t_Lattice::sites )
     .def_readwrite( "scale", &t_Lattice::scale )
     .def( "__str__",  &Ising_CE::details::print<t_Lattice> )
     .def( "fromXML",  &Ising_CE::details::fromXML<t_Lattice> );

   class_< t_CS >( "CS" )
     .def( init< t_CS >() )
     .def( "evaluate", &t_CS::evaluate )
     .def( "fromXML",  &Ising_CE::details::fromXML<t_CS> );

   class_< t_CH >( "ConvexHull" )
     .def( "__str__",  &t_CH::print )
     .def( "evaluate", &t_CH::evaluate )
     .def( "add",      &t_CH::add );

   class_< t_PiStructure >( "PiStructure" )
     .def( init< t_PiStructure >() )
     .def( init< types::t_int >() )
     .def( init< types::t_int, types::t_real >() )
     .def( "__str__",  &t_PiStructure::print )
     .def_readwrite( "index", &t_PiStructure::index )
     .def_readwrite( "x", &t_PiStructure::x );

   class_< t_Builder >( "Builder" )
     .def( init< t_Builder >() )
     .def( "fromXML",  &Ising_CE::details::fromXML<t_Builder> )
     .def( "build",    &Ising_CE::details::generateCEs,
           return_value_policy<manage_new_object>() );

   class_< t_Builder::t_Chemical >( "Chemical" )
     .def( init< t_Builder :: t_Chemical >() )
     .def( "evaluate", &t_Builder::t_Chemical::evaluate );

   class_< t_Builder::t_CS >( "ConstituentStrain" )
     .def( init< t_Builder :: t_CS >() )
     .def( "evaluate", &t_Builder::t_CS::evaluate );

   class_< CEFunc >( "CE", noinit )
     .def( "evaluate", &Ising_CE::details::CEFunc::evaluate );
     .def( "__len__",     &Ising_CE::details::CEFunc::size )
     .def( "__getitem__", &Ising_CE::details::CEFunc::getitem )
     .def( "__setitem__", &Ising_CE::details::CEFunc::setitem );
#endif

#include "atom.py.hpp"

#endif 

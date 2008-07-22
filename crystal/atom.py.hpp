//
//  Version: $Id$
//
using namespace Crystal;
using namespace CE;

typedef Crystal::Lattice     t_Lattice;
typedef Crystal::Structure   t_Structure;
typedef t_Structure::t_Atom   t_Atom;
typedef t_Lattice::t_Site     t_Site;


namespace details
{
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

  template< class T_HARMONIC >
  class CEFunc : public Builder<T_HARMONIC> :: t_VA_Functional 
  {
    typedef Builder<T_HARMONIC> t_Builder;
    typedef typename t_Builder :: t_VA_Functional t_Base;
    typedef typename t_Base :: t_Functional1 t_Chemical;
    typedef typename t_Base :: t_Functional2 t_CS;

    public:
      CEFunc() : t_Base() {} 
      CEFunc( t_Chemical* _chem, t_CS* _cs ) : t_Base( _chem, _cs ) {}
      ~CEFunc ()
      { 
        if( t_Base::functional1 ) delete t_Base::functional1;
        if( t_Base::functional2 ) delete t_Base::functional2;
      } 
      void assign( t_Structure &_str )
      {
        std::transform(
           _str.atoms.begin(), _str.atoms.end(), t_Base::begin(),
           boost::lambda::bind( &t_Structure::t_Atom::type, boost::lambda::_1 ) 
        );
      }
  };

  template< class T_TYPE > std::string print( const T_TYPE &_at );

    
  list symmetry_ops(  const t_Lattice& );



  t_Atom* AtomFromObject( list& _ob );
  t_Site* SiteFromObject( list& _ob );

  namespace XML
  {
    template< class T_TYPE > std::string nodename() { return "Unknown Type"; }
    template< class T_TYPE > void do_specialcode( T_TYPE &_type ) {}
    template< class T_TYPE > bool doloadcode( T_TYPE &_type, TiXmlElement *_parent )
      { return _type.Load( *_parent ); }
    template< class T_TYPE > TiXmlElement* findnode( TiXmlHandle &_doc )
      { return _doc.FirstChild("Job").FirstChild(nodename<T_TYPE>()).Element(); }
    template< class T_TYPE >
    void from(T_TYPE &_type, const std::string &_filename )
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

    template< class T_TYPE >
      void to( const T_TYPE &_type, const std::string &_filename )
      {
        
        TiXmlElement* parent = new TiXmlElement( nodename<T_TYPE>() ); 
        _type.print_xml( *parent );

        std::ifstream doesexist;
        doesexist.open(_filename.c_str(), std::ifstream::in);
        doesexist.close();
        TiXmlDocument doc;
        if(doesexist.fail())
        {
          doc.SetTabSize(1);
          doc.LinkEndChild( new TiXmlDeclaration("1.0", "", "") );
          TiXmlElement *node = new TiXmlElement("Job");
          node->LinkEndChild( parent );
          doc.LinkEndChild( node );
          doesexist.clear(std::ios::failbit);
        }
        else
        {
          doc.LoadFile( _filename.c_str());
          TiXmlHandle docHandle( &doc );
          TiXmlElement *child = docHandle.FirstChild("Job").Element();
          child->LinkEndChild( parent );
        }
        doc.SaveFile(_filename.c_str() );
      }
  }


  template<class T_HARMONIC>
  void generateCS( ConstituentStrain::Functional<T_HARMONIC> &_cs, t_Structure &_str )
  {
    __DOASSERT( _str.atoms.size() < 1, 
                "Cannot create constituent strain from "
                "structure with no atoms.\n" )
    _str.find_k_vectors();
    _cs << _str;
  }


  types::t_real toReal(std::string _str );
  std::string toType( types::t_real _r );

  template< class T_HARMONIC > 
  CEFunc<T_HARMONIC>* generateCEs( Builder<T_HARMONIC> &_builder, t_Structure &_str );

  template< class  T_FUNC > void assign( T_FUNC&, t_Structure&); 
  template< class  T_FUNC > void createCS( T_FUNC& _f, t_Structure& _s)
  { 
     _s.find_k_vectors();
    _f << _s; _f.resize( _s.atoms.size() ); 
  }  
  template< class  T_FUNC > typename T_FUNC::t_Container& get_vars( T_FUNC& _f )
  { return *_f.get_variables(); }


  template< class T_HARMONIC >
    void ExposeHarmonicRelated()
    {
      typedef T_HARMONIC t_Harmonic;
      typedef ConstituentStrain::Functional< t_Harmonic > t_CS;
      typedef Builder< t_Harmonic > t_Builder;
      typedef CEFunc< t_Harmonic > t_CEFunc;

      typename t_CS::t_Container* (t_CS::*varfunc)() const = &t_CS::get_variables;
      std::string name = t_Harmonic::type + "CS";
      class_< t_CS >( name.c_str() )
        .def( init< t_CS >() )
        .def( "evaluate", &t_CS::evaluate )
        .def( "assign",   &assign<t_CS> )
        .def( "define",   &createCS<t_CS> )
        .def( "vars",     varfunc,
              return_internal_reference<1>() )
        .def( "fromXML",  &XML::from< t_CS > );
     
      name = t_Harmonic::type + "Builder";
      class_< t_Builder >( name.c_str() )
        .def( init< t_Builder >() )
        .def( "fromXML",  &XML::from<t_Builder> )
        .def( "build",    &generateCEs<t_Harmonic>,
              return_value_policy<manage_new_object>() );
      
      name = t_Harmonic::type + "CE";
      class_< t_CEFunc >( name.c_str(), no_init )
        .def( "evaluate", &t_CEFunc::evaluate )
        .def( "assign",   &t_CEFunc::assign )
        .def( "vars",     &t_CEFunc::get_variables,
              return_internal_reference<1>() )
        .def( "chemical", &t_CEFunc::get_functional1,
              return_internal_reference<1>() )
        .def( "CS", &t_CEFunc::get_functional2,
              return_internal_reference<1>() );
    }

# include "atom.py.impl.hpp"
}



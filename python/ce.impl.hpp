//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/python.hpp>

#include <ce/functional_builder.h>
#include <ce/constituent_strain.h>
#include <ce/harmonic.h>
#include <crystal/lattice.h>
#include <crystal/structure.h>
#include <crystal/atom.h>

#include "xml.hpp"


namespace LaDa
{
  namespace Python
  {
    namespace XML
    {
      typedef CE::Builder< CE::ConstituentStrain::Harmonic::Cubic > t_CubicBuilder;
      typedef t_CubicBuilder :: t_CS t_CubicCS;
      typedef CE::Builder< CE::ConstituentStrain::Harmonic::Tetragonal > t_TetraBuilder;
      typedef t_TetraBuilder :: t_CS t_TetraCS;

      template<> std::string nodename<t_CubicCS>();
      template<> std::string nodename<t_TetraCS>();
      
      template<> void do_specialcode< t_CubicBuilder >( t_CubicBuilder &_type );
      template<> void do_specialcode< t_TetraBuilder >( t_TetraBuilder &_type );
      
      template<> bool doloadcode<t_CubicCS>( t_CubicCS &_type, TiXmlElement *_parent );
      template<> bool doloadcode<t_TetraCS>( t_TetraCS &_type, TiXmlElement *_parent );
        
      template<> TiXmlElement *findnode<t_CubicCS>( TiXmlHandle &_doc );
      template<> TiXmlElement *findnode<t_TetraCS>( TiXmlHandle &_doc );
      template<> TiXmlElement *findnode<t_CubicBuilder>( TiXmlHandle &_doc );
      template<> TiXmlElement *findnode<t_TetraBuilder>( TiXmlHandle &_doc );
    }

    void export_ce();

    namespace details
    {
   
      template< class T_HARMONIC >
      class CEFunc : public CE::Builder<T_HARMONIC> :: t_VA_Functional 
      {
        typedef CE::Builder<T_HARMONIC> t_Builder;
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
          void assign( Crystal::Structure &_str )
          {
            std::transform(
               _str.atoms.begin(), _str.atoms.end(), t_Base::begin(),
               boost::lambda::bind( &Crystal::Structure::t_Atom::type, boost::lambda::_1 ) 
            );
          }
      };

      template<class T_HARMONIC>
      void generateCS( CE::ConstituentStrain::Functional<T_HARMONIC> &_cs,
                       Crystal::Structure &_str )
      {
        __DOASSERT( _str.atoms.size() < 1, 
                    "Cannot create constituent strain from "
                    "structure with no atoms.\n" )
        _str.find_k_vectors();
        _cs << _str;
      }
   
      template< class T_HARMONIC > 
      CEFunc<T_HARMONIC>* generateCEs( CE::Builder<T_HARMONIC> &_builder,
                                       Crystal::Structure &_str );
   
      template< class  T_FUNC > void assign( T_FUNC&, Crystal::Structure&); 
      template< class  T_FUNC > void createCS( T_FUNC& _f, Crystal::Structure& _s)
      { 
         _s.find_k_vectors();
         foreach( Crystal::Structure::t_kAtom &kat, _s.k_vecs )
        _f << _s; _f.resize( _s.atoms.size() ); 
      }  
      template< class  T_FUNC > typename T_FUNC::t_Container& get_vars( T_FUNC& _f )
      { return *_f.get_variables(); }
   
      template< class T_HARMONIC > 
      CEFunc<T_HARMONIC>* generateCEs( CE::Builder<T_HARMONIC> &_builder, Crystal::Structure &_str )
      {
        typedef CE::Builder<T_HARMONIC> t_Builder;
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
      void assign( T_FUNC &_func, Crystal::Structure &_str )
      {
        std::transform(
           _str.atoms.begin(), _str.atoms.end(), _func.begin(),
           boost::lambda::bind( &Crystal::Structure::t_Atom::type, boost::lambda::_1 ) 
        );
      }

      template< class T_HARMONIC >
        void ExposeHarmonicRelated()
        {
          using namespace boost::python;
          typedef T_HARMONIC t_Harmonic;
          typedef CE::ConstituentStrain::Functional< t_Harmonic > t_CS;
          typedef CE::Builder< t_Harmonic > t_Builder;
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

    }
  }
} // namespace LaDa

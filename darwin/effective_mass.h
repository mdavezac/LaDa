//
//  Version: $Id$
//
#ifndef _DARWIN_EFFECTIVE_MASS_STUBS_H_
#define _DARWIN_EFFECTIVE_MASS_STUBS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/serialization/string.hpp>

#include <tinyxml/tinyxml.h>

namespace LaDa
{
  namespace GA
  {
    namespace Keepers
    {
#     define DECLARE_KEEPER( type, classname, name ) \
      struct classname \
      { \
        friend class boost::serialization::access; \
        type name; \
        classname() : name(0) {} \
        classname( const classname& _c ) : name(_c.name) {}\
        bool Load( const TiXmlElement &_node )\
        {\
          if ( not _node.Attribute(#name ) )\
          {\
            std::cerr << "Could not Load GA::Keepers::" << #classname << ".\n";\
            return false;\
          }\
          name = boost::lexical_cast<type>( _node.Attribute(#name) );\
          return true;\
        }\
        bool Save( TiXmlElement &_node ) const\
        {\
          const std::string value = boost::lexical_cast<std::string>( name );\
          _node.SetAttribute( #name, value );\
          return true;\
        }\
        private:\
          template<class Archive> void serialize(Archive & _ar, const unsigned int _version)\
            { _ar & name; }\
      };

      //! Object stub which keeps track of electron effective mass.
      DECLARE_KEEPER( types::t_real, eMass, emass );
      //! Object stub which keeps track of hole effective mass.
      DECLARE_KEEPER( types::t_real, hMass, hmass );

#     undef DECLARE_KEEPER
    } // namespace Keepers
  } // namespace GA

} // namespace LaDA

#endif // _PESCAN_H_

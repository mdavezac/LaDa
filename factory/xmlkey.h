//
//  Version: $Id: factory.h 850 2008-11-13 00:53:54Z davezac $
//
#ifndef _LADA_FACTORY_XMLKEY_H_
#define _LADA_FACTORY_XMLKEY_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<string>
#include<sstream>
#include<utility>

namespace LaDa
{
  namespace Factory
  {
    //! Key for XML factory.
    struct XmlKey 
    {
      //! Constructor.
      XmlKey   ( const std::string& _a, const std::string& _b = "", const std::string _c = "" )
             : node( _a ), attribute(_b), value( _c )  {}
      //! Copy Constructor.
      XmlKey   ( const XmlKey& _c )
             : node( _c.node ), attribute(_c.attribute ), value( _c.value )  {}

      //! Node name.
      const std::string node;
      //! Attribute name.
      const std::string attribute;
      //! Value name.
      const std::string value;

      //! Converts to a string.
      operator const std::string() const;
    };

    //! Compares two xml keys.
    inline bool operator < ( const XmlKey &_a, const XmlKey &_b )
    {
      if( _a.node != _b.node ) return _a.node < _b.node; 
      if( _a.attribute != _b.attribute ) return _a.attribute < _b.attribute; 
      return _a.value < _b.value; 
    }
    //! Compares two xml keys.
    inline bool operator > ( const XmlKey &_a, const XmlKey &_b ) 
    {
      if( _a.node != _b.node ) return _a.node > _b.node; 
      if( _a.attribute != _b.attribute ) return _a.attribute > _b.attribute; 
      return _a.value > _b.value; 
    }
    //! Compares two xml keys.
    inline bool operator == ( const XmlKey &_a, const XmlKey &_b ) 
    {
      return     _a.node == _b.node
             and _a.attribute == _b.attribute
             and _a.value == _b.value;
    }
    //! Dumps a key to a string.
    XmlKey::operator const std::string() const
    {
      std::ostringstream ss;
      if( attribute.empty() and value.empty() )
        return "<" + node + "/>";
      ss << "<" << ( node.empty() ? "??": node )
         << " " << ( attribute.empty() ? "??": attribute )
         << "=\"" << ( value.empty() ? "??" : value ) << "\" />";

      return ss.str();
    }
    inline std::ostream& operator<<( std::ostream& _stream, const XmlKey& _key )
      { return _stream << (const std::string) _key; }

  } // namespace Factory.
} // namespace LaDa.

#endif

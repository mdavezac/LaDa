//
//  Version: $Id: factory.h 860 2008-11-17 18:37:10Z davezac $
//
#ifndef _LADA_FACTORY_VISIT_XML_H_
#define _LADA_FACTORY_VISIT_XML_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <tinyxml/tinyxml.h>
#include <opt/tinyxml.h>

namespace LaDa
{
  //! Holds factory related objects.
  namespace Factory
  {
    //! Visits all children of a node and applies a factory.
    template< class T_FACTORY, class T_ARG > 
      void visit_xml( const T_FACTORY &_factory, const TiXmlElement &_node,
                      T_ARG& _arg, size_t _reentrancy = 1  )
      {
        typedef T_FACTORY t_Factory;
        typedef typename t_Factory::t_Key t_Key;

        const TiXmlElement *child = _node.FirstChildElement();
        for(; child; child = child->NextSiblingElement() )
        {
          const std::string name = child->Value();
          // Tries without attributes first.
          _factory( t_Key( name ), _arg, *child );

          // Then tries attributes.
          opt::const_AttributeIterator i_att( *child );
          opt::const_AttributeIterator i_att_end;
          for(; i_att != i_att_end; ++i_att )
          {
            // attribute alone. no node name.
            _factory( t_Key("", i_att->first ), _arg, i_att->second );

            // then attribute + name
            _factory( t_Key(name, i_att->first ), _arg, i_att->second );
            
            // then attribute + name + value
            _factory( t_Key(name, i_att->first, i_att->second ), _arg );
          } // loop over attributes.

          if( _reentrancy > 1 ) visit_xml( _factory, *child, _arg, _reentrancy - 1 );
          if( _reentrancy == 0 ) visit_xml( _factory, *child, _arg, _reentrancy );
        } // loop over elements.
      }

  } // namespace factory.
}

#endif 

#ifndef __PYTHONLADA_XML_HPP_
#define __PYTHONLADA_XML_HPP_

#include "LaDaConfig.h"

#include<fstream>

#include<tinyxml/tinyxml.h>

#include<opt/types.h>
#include<opt/debug.h>

namespace LaDa
{
  namespace Python
  {
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
        do_specialcode<T_TYPE>( _type );
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
  }
} // namespace LaDa

#endif // __PYTHONLADA_XML_HPP_

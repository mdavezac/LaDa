//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/filesystem/operations.hpp>

#include <tinyxml/tinyxml.h>

#include "save_every.h"


namespace LaDa
{
  namespace GA
  {
    namespace CheckPoint
    {
      void SaveGAState :: operator()(bool) const
      {
        namespace bfs = boost::filesystem;
        __TRYBEGIN
          TiXmlDocument doc;
          TiXmlElement *node;
          bool reopen = bfs::exists( filename_ );
          if ( reopen )
          {
            doc.LoadFile(filename_.string());
            TiXmlHandle docHandle( &doc );
            node = docHandle.FirstChild("Job").Element();
            if ( not node ) reopen = false;
            {
              TiXmlNode* child = docHandle.FirstChild("Job")
                                          .FirstChild("Restart").Node();
              for(; child; child = docHandle.FirstChild("Job").FirstChild("Restart").Node() )
                node->RemoveChild(child);
              node = new TiXmlElement("Restart");
              __DOASSERT( not node, "Memory allocation error.\n" )
              child = docHandle.FirstChild("Job").Element();
              child->LinkEndChild( node );
            }
          }
          if( not reopen )
          {
            doc.SetTabSize(1);
            doc.LinkEndChild( new TiXmlDeclaration("1.0", "", "") );
            TiXmlElement *jobnode = new TiXmlElement("Job");
            __DOASSERT( not jobnode, "Memory allocation error.\n" )
            doc.LinkEndChild( jobnode );
            node = new TiXmlElement("Restart");
            __DOASSERT( not node, "Memory allocation error.\n" )
            jobnode->LinkEndChild( node );
          }

          (*signals_)( *node );
          __DOASSERT( not doc.SaveFile(filename_.string() ), "" )
        __ENDGROUP__
        catch(...)
        {
          std::cerr << "Error while trying to saved GA state.\n";
          Print::out << "Error while trying to save GA state.\n";
          return;
        }
        Print::out << "Save GA state in " << filename_ << Print::endl;
      } 
    } // namespace CheckPoint
  } // namespace GA
} // namespace LaDa

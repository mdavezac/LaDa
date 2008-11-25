//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<boost/filesystem/operations.hpp>
#include<boost/mpi/collectives.hpp>
#include<boost/bind.hpp>
#include<boost/ref.hpp>
#include<iterator>
#include<fstream>
#include<sstream>

#include <mpi/macros.h>
#include <mpi/mpi_object.h>
#include <opt/debug.h>
#include "tinyxml.h"

namespace LaDa
{
  namespace opt
  {
    //! \brief Returns the node \a _name.
    //! \details Looks first to \a _element, then its childrent, then its
    //!          next siblings.
    //! \todo Look amongst all siblings.
    const TiXmlElement* find_node( const TiXmlElement &_element,
                                   const std::string& _name )
    {
      const TiXmlElement *parent;
    
      // Find first XML "Structure" node (may be _element from start, or a child of _element)
      std::string str = _element.Value();
      if ( str.compare(_name) == 0 ) return &_element;
      parent =  _element.FirstChildElement(_name);
      if( parent ) return parent;

      return _element.NextSiblingElement(_name);
    }


    //! \brief Returns the node \<Functional type=\a _name\>.
    //! \details Looks first to \a _element, then its childrent, then its
    //!          next siblings.
    const TiXmlElement* find_functional_node ( const TiXmlElement &_element,
                                               const std::string &_name )
    {
      const TiXmlElement *parent;
      std::string str;

      // This whole section tries to find a <Functional type="vff"> tag
      // in _element or its child
      str = _element.Value();
      if ( str.compare("Functional" ) != 0 )
        parent = _element.FirstChildElement("Functional");
      else parent = &_element;
      
      while (parent)
      {
        str = "";
        if ( parent->Attribute( "type" )  )
          str = parent->Attribute("type");
        if ( str.compare(_name) == 0 )
          break;
        parent = parent->NextSiblingElement("Functional");
      }
      if ( parent ) return parent;
      
      std::cerr << "Could not find an <Functional type=\""
                << _name << "\"> tag in input file" 
                << std::endl;
      return NULL;
    }  // Functional :: find_node


    void read_file( const boost::filesystem::path &_input, std::string& _result )
    {
      namespace bfs = boost::filesystem;
      __ROOTCODE
      ( 
        (*::LaDa::mpi::main), 
         __DOASSERT( not bfs::exists( _input ), _input << " does not exist.\n" )
         __DOASSERT( not ( bfs::is_regular( _input ) or bfs::is_symlink( _input ) ),
                     _input << " is a not a valid file.\n" );
        
         std::ifstream file( _input.string().c_str(), std::ios_base::in );
         std::istream_iterator< std::string > i_file( file );
         std::istream_iterator< std::string > i_file_end;
         for(; i_file != i_file_end; ++i_file ) _result.append( *i_file );
         file.close();
      )
      __MPICODE( boost::mpi::broadcast( *LaDa::mpi::main, _result, 0 ); )
    }
    void read_xmlfile( const boost::filesystem::path &_input, std::string& _result )
    {
      __ROOTCODE
      ( 
        (*::LaDa::mpi::main), 
        TiXmlDocument doc( _input.string() ); 
        TiXmlHandle docHandle( &doc ); 
        __DOASSERT( not doc.LoadFile(), 
                       doc.ErrorDesc() << "\n" 
                    << "Could not load input file " << _input 
                    << " in  Darwin<T_EVALUATOR>.\nAborting.\n" ) 
        std::ostringstream stream;
        const TiXmlElement *parent = doc.RootElement();
        stream << *parent;
        _result = stream.str();
      )
      __MPICODE( boost::mpi::broadcast( *LaDa::mpi::main, _result, 0 ); )
    }

  }
} // namespace LaDa

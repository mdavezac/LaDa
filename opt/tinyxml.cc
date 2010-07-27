#include "LaDaConfig.h"

#include<iterator>
#include<fstream>
#include<sstream>

#include<boost/filesystem/path.hpp>
#include<boost/filesystem/operations.hpp>
#include<boost/bind.hpp>
#include<boost/ref.hpp>
#ifdef LADA_MPI
# include<boost/mpi/collectives.hpp>
#endif

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
                                   const std::string& _name,
                                   const std::string& _attribute,
                                   const std::string& _value )
    {
      const TiXmlElement *parent;
    
      // Find first XML "Structure" node (may be _element from start, or a child of _element)
      std::string str = _element.Value();
      if ( str.compare(_name) == 0 )
      {
        if( _attribute == "" ) return &_element;
        if( _element.Attribute(_attribute) )
        {
          if( _value == "" ) return &_element;
          const std::string value = *_element.Attribute( _attribute );
          if( _value == value ) return &_element;
        }
      }
      parent =  _element.FirstChildElement(_name);
      if( parent and _attribute == "" and _value == ""  ) return parent;
      for(; parent; parent = parent->NextSiblingElement( _name ) )
      {
        if( not parent->Attribute( _attribute ) ) continue;
        if( _value == "" ) return parent;
        const std::string value = *parent->Attribute( _attribute );
        if( _value == value ) return parent;
      }

      parent = _element.NextSiblingElement(_name);
      if( parent and _attribute == "" and _value == ""  ) return parent;
      for(; parent; parent = parent->NextSiblingElement( _name ) )
      {
        if( not parent->Attribute( _attribute ) ) continue;
        if( _value == "" ) return parent;
        const std::string value = *parent->Attribute( _attribute );
        if( _value == value ) return parent;
      }

      return NULL;
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
      LADA_MPI_CODE( boost::mpi::communicator world; )
      LADA_ROOT
      ( 
        world, 
         LADA_DO_NASSERT( not bfs::exists( _input ), _input << " does not exist.\n" )
         LADA_DO_NASSERT( not ( bfs::is_regular( _input ) or bfs::is_symlink( _input ) ),
                     _input << " is a not a valid file.\n" );
        
         std::ifstream file( _input.string().c_str(), std::ios_base::in );
         std::istream_iterator< std::string > i_file( file );
         std::istream_iterator< std::string > i_file_end;
         for(; i_file != i_file_end; ++i_file ) _result.append( *i_file );
         file.close();
      )
      LADA_MPI_CODE( boost::mpi::broadcast( world, _result, 0 ); )
    }
    void read_xmlfile( const boost::filesystem::path &_input, std::string& _result )
    {
      LADA_MPI_CODE( boost::mpi::communicator world; )
      LADA_ROOT
      ( 
        world, 
        TiXmlDocument doc( _input.string() ); 
        TiXmlHandle docHandle( &doc ); 
        LADA_DO_NASSERT( not doc.LoadFile(), 
                       doc.ErrorDesc() << "\n" 
                    << "Could not load input file " << _input << ".\n" )
        std::ostringstream stream;
        const TiXmlElement *parent = doc.RootElement();
        stream << *parent;
        _result = stream.str();
      )
      LADA_MPI_CODE( boost::mpi::broadcast( world, _result, 0 ); )
    }
    void read_xmlfile( const boost::filesystem::path &_input, TiXmlDocument& _doc )
    {
      std::string string;
      read_xmlfile( _input, string );
      _doc.Parse( string.c_str() );
    }

  }
} // namespace LaDa

//
//  Version: $Id$
//
#ifndef _OPT_TINYXML_H_
#define _OPT_TINYXML_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/filesystem/path.hpp>
#include <string>
#include <tinyxml/tinyxml.h>
#include <opt/path.h>
#include "debug.h"

namespace LaDa
{
  namespace opt
  {
    //! \brief Returns the node \a _name.
    //! \details Looks first to \a _element, then its children, then its
    //!          next siblings. 
    //! \todo Look amongst all siblings.
    const TiXmlElement* find_node( const TiXmlElement &_element,
                                   const std::string& _name,
                                   const std::string& _attribute = "",
                                   const std::string& _value = "" );


    //! \brief Returns the node \<Functional type=\a _name\>.
    //! \details Looks first to \a _element, then its childrent, then its
    //!          next siblings.
    const TiXmlElement* find_functional_node ( const TiXmlElement &_element,
                                               const std::string &_name );

    //! \brief reads a file and dumps it in a string.
    //! \details Works in mpi.
    void read_file( const boost::filesystem::path &_input, std::string& _result );
    //! \brief reads an file and dumps it in a string.
    //! \details Works in mpi.
    void read_xmlfile( const boost::filesystem::path &_input, std::string& _result );
    //! \brief reads an file and has a TiXmlDocument parse it.
    void read_xmlfile( const boost::filesystem::path &_input, TiXmlDocument& _doc );
      
    //! \brief reads a tag from xml, with filename redirection.
    template< class T_FUNCTIONAL >
      void read_tag( T_FUNCTIONAL &_functional, 
                     const TiXmlElement &_node, 
                     const std::string &_name,
                     const std::string& _attribute = "",
                     const std::string& _value = ""  );
    //! \brief reads a tag from input file.
    template< class T_FUNCTIONAL >
      void read_tag( T_FUNCTIONAL &_functional,
                     const boost::filesystem::path &_path,
                     const std::string &_name,
                     const std::string& _attribute = "",
                     const std::string& _value = ""  );

    //! Constant attribute iterator.
    class const_AttributeIterator 
    {
      public:
        //! return type on dereference.
        typedef std::pair< std::string, std::string > value_type;

        //! Constructor.
        const_AttributeIterator() : att_(NULL), node_(NULL) {}
        //! Constructor.
        const_AttributeIterator   ( const TiXmlElement& _node )
           { set_node( _node ); }
        //! Constructor.
        const_AttributeIterator   ( const const_AttributeIterator &_it )
                                : att_(_it.att_), node_( _it.node_ ) {}
        //! Dereference.
        const value_type operator*() const
        { 
          __ASSERT( not ( node_ or att_ ), "iterator points to nothing." )
          return value_type( att_->Name(), att_->Value() );
        }
        //! Dereference.
        const value_type* operator->() const
        { 
          __ASSERT( not ( node_ or att_ ), "iterator points to nothing." )
          buffer.first = att_->Name(); 
          buffer.second = att_->Value(); 
          return &buffer;
        }

        //! Pre-increment.
        const_AttributeIterator operator++() 
        { 
          __ASSERT( not ( node_ or att_ ), "iterator points to nothing." )
          att_ = att_->Next();
          return *this;
        }
        //! Post-increment.
        const_AttributeIterator operator++( int ) 
        { 
          __ASSERT( not ( node_ or att_ ), "iterator points to nothing." )
          const_AttributeIterator old( *this );
          att_ = att_->Next();
          return old;
        }
        //! Equality.
        bool operator == ( const const_AttributeIterator& _a ) const
        {
          return    ( node_ == _a.node_ and att_ == _a.att_ )
                 or ( _a.att_ == NULL and att_ == NULL );
        }
        //! Non-equality.
        bool operator != ( const const_AttributeIterator& _a ) const
        {
          return     ( node_ != _a.node_ or att_ != _a.att_ )
                 and ( _a.att_ or att_ );
        }


        //! Sets the node.
        void set_node( const TiXmlElement &_node )
          { node_ = &_node; att_ = node_->FirstAttribute(); }
        //! Gets the node.
        const TiXmlElement& get_node() const { return *node_; }

        const TiXmlAttribute* as_tinyxml() const { return att_; }

      protected:
        //! Pointer to current position.
        const TiXmlAttribute *att_; 
        //! Pointer to current node.
        const TiXmlElement *node_; 
        //! An internal buffer.
        mutable value_type buffer;
    };

    //! Explores all the attributes of \a _node.
    template< class T_FACTORY >
      void explore_attributes( T_FACTORY& _factory, const TiXmlElement &_node )
      {
        opt::const_AttributeIterator i_att( _node );
        const opt::const_AttributeIterator i_att_end;
        for(; i_att != i_att_end; ++i_att )
          _factory( i_att->first, i_att->second );
      }

//   template< class T_FUNCTIONAL >
//     void read_functional( T_FUNCTIONAL &_functional,
//                           const TiXmlElement &_node,
//                           const std::string& _name  )
//     {
//       const TiXmlElement* node = find_functional_node( _node, _name );
//       __DOASSERT( not node, "Could not read functional " + _name + " from input.\n" )
//       if( node->Attribute("filename") )
//         read_functional( _functional,
//                          boost::filesystem::path
//                          (
//                            opt::expand_path(node->Attribute("filename"))
//                          ),  _name );
//       _functional.Load( *node );
//     }
//
//   template< class T_FUNCTIONAL >
//     void read_functional( T_FUNCTIONAL &_functional,
//                           const boost::filesystem::path &_path,
//                           const std::string& _name  )
//     {
//       TiXmlDocument doc;
//       TiXmlHandle docHandle( &doc ); 
//       read_xmlfile( _path, doc );
//       const TiXmlElement* child = docHandle.FirstChild( "Job" ).Element();
//       __DOASSERT( not child, "Could not find root node \"Job\".\n" )
//       read_functional( _functional, *child, _name );
//     }

    template< class T_OBJECT >
      void read_tag( T_OBJECT &_object,
                     const TiXmlElement &_node,
                     const std::string& _name,
                     const std::string& _attribute,
                     const std::string& _value )
      {
        const TiXmlElement* node = find_node( _node, _name, _attribute, _value );
        const std::string tagname
        (
            "<" + _name 
          + ( 
              _attribute == ""  ?
              "":
              " " + _attribute + ( _value == ""  ? "=\"??\"": "=\"" + _value + "\"" )
            ) 
          + "/>"
        );
        __DOASSERT( not node, "Could not read tag " + tagname + " from input.\n" )
        if( node->Attribute("filename") )
        {
          read_tag( _object,
                    boost::filesystem::path
                    ( 
                      opt::expand_path(node->Attribute("filename"))
                    ),  _name, _attribute, _value );
          return;
        }
        _object.Load( *node );
      }

    template< class T_OBJECT >
      void read_tag( T_OBJECT &_object,
                     const boost::filesystem::path &_path,
                     const std::string& _name,
                     const std::string& _attribute,
                     const std::string& _value )
      {
        TiXmlDocument doc;
        TiXmlHandle docHandle( &doc ); 
        read_xmlfile( _path, doc );
        const TiXmlElement* child = docHandle.FirstChild( "Job" ).Element();
        __DOASSERT( not child, "Could not find root node \"Job\".\n" )
        read_tag( _object, *child, _name, _attribute, _value );
      }
  }
} // namespace LaDa
#endif

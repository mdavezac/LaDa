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
#include "debug.h"

namespace LaDa
{
  namespace opt
  {
    //! \brief Returns the node \a _name.
    //! \details Looks first to \a _element, then its childrent, then its
    //!          next siblings.
    //! \todo Look amongst all siblings.
    const TiXmlElement * find_node( const TiXmlElement &_element,
                                    const std::string& _name );

    //! \brief Returns the node \<Functional type=\a _name\>.
    //! \details Looks first to \a _element, then its childrent, then its
    //!          next siblings.
    const TiXmlElement* find_functional_node ( const TiXmlElement &_element,
                                               const std::string &_name );

    //! \brief reads a file and dumps it in a string.
    //! \details Works in mpi.
    void read_file( boost::filesystem::path &_input, std::string& _result );
    //! \brief reads an file and dumps it in a string.
    //! \details Works in mpi.
    void read_xmlfile( boost::filesystem::path &_input, std::string& _result );

    //! Attribute iterator.
    class const_AttributeIterator 
    {
      public:
        //! return type on dereference.
        typedef const std::pair< std::string, std::string > value_type;

        //! Constructor.
        const_AttributeIterator() : att_(NULL), node_(NULL) {}
        //! Constructor.
        const_AttributeIterator   ( const TiXmlElement& _node )
           { set_node( _node ); }
        //! Constructor.
        const_AttributeIterator   ( const const_AttributeIterator &_it )
                                : att_(_it.att_), node_( _it.node_ ) {}
        //! Dereference.
        value_type operator*()
        { 
          __ASSERT( not ( node_ or att_ ), "iterator points to nothing." )
          return value_type( att_->Name(), att_->Value() );
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
        bool operator == ( const const_AttributeIterator& _a )
        {
          return    ( node_ == _a.node_ and att_ == _a.att_ )
                 or ( _a.att_ == NULL and att_ == NULL );
        }

        //! Sets the node.
        void set_node( const TiXmlElement &_node )
          { node_ = &_node; att_ = node_->FirstAttribute(); }
        //! Gets the node.
        const TiXmlElement& get_node() const { return *node_; }

      protected:
        //! Pointer to current position.
        const TiXmlAttribute *att_; 
        //! Pointer to current node.
        const TiXmlElement *node_; 
    };

  }
} // namespace LaDa
#endif

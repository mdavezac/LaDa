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

  }
} // namespace LaDa
#endif

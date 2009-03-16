//
//  Version: $Id$
//

#ifndef __PYTHONLADA_MISC_HPP_
#define __PYTHONLADA_MISC_HPP_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <sstream>
#include <string>

namespace LaDa
{
  namespace Python
  {
    template< class T_TYPE >
    std::string print( const T_TYPE &_at )
    { 
      std::ostringstream sstr;
      _at.print_out( sstr );
      return sstr.str();
    }
    template< class T_TYPE >
    std::string tostream( const T_TYPE &_at )
    { 
      std::ostringstream sstr;
      sstr << _at;
      return sstr.str();
    }
    template< class T_TYPE, template< class > class T_CONTAINER > 
      void exposeContainer( std::string &_name )
      {
        namespace bp = boost::python;
        bp::class_< T_CONTAINER<T_TYPE> >( _name.c_str() )
          .def(bp::vector_indexing_suite< T_CONTAINER< T_TYPE > >());
      }

  }
} // namespace LaDa
#endif // __PYTHONLADA_MISC_HPP_

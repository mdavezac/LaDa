//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/filesystem/path.hpp>
#include <boost/python/def.hpp>
#include <boost/python/errors.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/manage_new_object.hpp>

#include "../read_structure.h"

namespace LaDa
{
  namespace Python
  {
    namespace details
    {
      Crystal::Structure* read_structure( const std::string &_path )
      {
        Crystal::Structure *result = new Crystal::Structure();
        try
        { 
          Crystal::read_structure( *result, _path ); 
          return result;
        }
        catch( std::exception &_e )
        {
          delete result;
          PyErr_SetString( PyExc_RuntimeError, _e.what() );
          boost::python::throw_error_already_set();
          return NULL;
        }
        catch( ... )
        {
          delete result;
          PyErr_SetString( PyExc_RuntimeError,
                           ("Could not read MBCE-type structure from " + _path + ".\n").c_str()  );
          boost::python::throw_error_already_set();
          return NULL;
        }
      }
    }

    void expose_read_structure()
    {
      boost::python::def
      (
        "read_structure", 
        &details::read_structure,
        boost::python::return_value_policy< boost::python::manage_new_object >(),
        "Tries to read a file as a MBCE structure file. Returns a structure on success" 
      );
    }

  }
} // namespace LaDa

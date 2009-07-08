//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <sstream>

#include <boost/filesystem/path.hpp>
#include <boost/python/def.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/manage_new_object.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/str.hpp>

#include "../read_structure.h"
#include <python/debug.hpp>

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
          LADA_PYTHON_ERROR( PyExc_RuntimeError, _e.what() );
          boost::python::throw_error_already_set();
          return NULL;
        }
        catch( ... )
        {
          delete result;
          LADA_PYTHON_ERROR
          ( 
            PyExc_RuntimeError,
            ("Could not read MBCE-type structure from " + _path + ".\n")
          );
          boost::python::throw_error_already_set();
          return NULL;
        }
      }

      void read_pifile_structure( boost::python::str& _str, Crystal::Structure &_structure )
      {
        try
        {
          const std::string string = boost::python::extract<std::string>( _str );
          if( string.empty() )
          {
            LADA_PYTHON_ERROR
            ( 
              PyExc_RuntimeError,
              ( "Could not read string.\n" )
            );
            boost::python::throw_error_already_set();
          }
          std::istringstream sstring( string );
          Crystal::read_pifile_structure( sstring, _structure );
        }
        catch( std::exception &_e )
        {
          LADA_PYTHON_ERROR( PyExc_RuntimeError, _e.what() );
          boost::python::throw_error_already_set();
        }
        catch( ... )
        {
          LADA_PYTHON_ERROR
          ( 
            PyExc_RuntimeError,
            "Could not read structure from pi-file format.\n"
          );
          boost::python::throw_error_already_set();
        }
        return;
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
      boost::python::def
      ( 
        "read_pifile_structure",
        &details::read_pifile_structure,
        ( boost::python::arg("input"), boost::python::arg("structure") ),
        "Reads a structure from pi-file type input.\n" 
      );

    }

  }
} // namespace LaDa

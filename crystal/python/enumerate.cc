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

#include "../enumerate.h"
#include <python/debug.hpp>

namespace LaDa
{
  namespace Python
  {
    namespace details
    {
      void enumerate_pifile( const std::string &_filename,
                             const boost::python::object &_callable )
      {
        try
        {
          Crystal :: Structure structure;
          std::ifstream file( _filename.c_str(), std::ifstream::in );
          do
          {
            if( not Crystal :: read_pifile_structure( file, structure ) ) continue;
            _callable( structure );
            foreach( Crystal::Structure::t_Atom &atom, structure.atoms )
              atom.type = Fuzzy::gt( atom.type, 0e0 ) ? -1e0: 1e0;
            structure.name = "-" + structure.name;
            _callable( structure );
          }
          while( not file.eof() );
        }
        catch( std::exception &_e )
        {
          LADA_PYTHON_ERROR
          ( 
            PyExc_RuntimeError,
            ( "Error while enumerating PI-file " + _filename + ":\n" + _e.what() )
          );
          boost::python::throw_error_already_set();
        }
        catch( ... )
        {
          LADA_PYTHON_ERROR
          ( 
            PyExc_RuntimeError,
            ( "Error while enumerating PI-file " + _filename + "\n" )
          );
          boost::python::throw_error_already_set();
        }
      }
    } // namespace details

    void expose_enumerate()
    {

      boost::python::def
      ( 
        "enum_pifile",
        &details::enumerate_pifile,
        ( boost::python::arg("filename"), boost::python::arg("callable") ),
        "Reads a pi-file and calls a "
        "\"callable( const LaDa::Crystal::Structure )->types::t_real\" for each structure."
      );
    }

  }
} // namespace LaDa

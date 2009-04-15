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
#include <boost/python/extract.hpp>

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
          PyErr_SetString
          ( 
            PyExc_RuntimeError,
            ("Could not read MBCE-type structure from " + _path + ".\n").c_str()  
          );
          boost::python::throw_error_already_set();
          return NULL;
        }
      }
      class EnumOp
      {
        public:
          EnumOp( const boost::python::object &_o ) : object_(_o) {}
          types::t_real operator()( const Crystal::Structure &_structure )
          {  
            try
            {
              const boost::python::object result( object_( _structure ) );
              const types::t_real r = boost::python::extract< types::t_real >( result );
              return r;
            }
            catch( ... )
            {
              PyErr_SetString(PyExc_RuntimeError,
                              "Error calling python from C, or extracting result.\n" );
              boost::python::throw_error_already_set();
            }
          }
        private:
          const boost::python::object &object_;
      };
      void enumerate_pifile( const std::string &_filename,
                             const boost::python::object &_callable )
      {
        try
        {
          EnumOp op( _callable );
          Crystal::enumerate_pifile( _filename, op );
        }
        catch( std::exception &_e )
        {
          PyErr_SetString( PyExc_RuntimeError, _e.what() );
          boost::python::throw_error_already_set();
        }
        catch( ... )
        {
          PyErr_SetString
          ( 
            PyExc_RuntimeError,
            ( "Error while enumerating PI-file " + _filename + "\n" ).c_str() 
          );
          boost::python::throw_error_already_set();
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

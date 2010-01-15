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

      void enumerate_gusfile( const std::string &_file, 
                              const boost::python::object &_callable )
      {
        try
        {
          Crystal :: Structure structure;
          __DOASSERT( structure.lattice == NULL, "Lattice not set in structure.\n" );
          std::ifstream file( _file.c_str(), std::ifstream::in );
          std::string line;
          do
          {
            getline(file, line);
            if( line.find("#tot") != std::string::npos ) break;
          } while( (not file.eof()) and file.good() );
         
          if( file.bad() ) return;
          do
          {
            getline(file, line);
            line = Print::StripEdges(line);
            if( line[0] == '#' ) continue;
            if( line.size() == 0 ) continue;
            std::istringstream input( line );
            structure.cell.zero();
            types::t_int dummy;
            input >> structure.name;
            for( size_t i(0); i < 6; ++i ) input >> dummy;
            input >> structure.cell(0,0) 
                  >> structure.cell(1,0) >> structure.cell(1,1)
                  >> structure.cell(2,0) >> structure.cell(2,1) >> structure.cell(2,2);
            for( size_t i(0); i < 9; ++i ) input >> dummy;
            input >> line;
            structure.cell = structure.lattice->cell * structure.cell;
            structure.atoms.clear();
            fill_structure( structure );
            __DOASSERT( line.size() != structure.atoms.size(),
                        "labels and structure have different sizes.\n" )
            Crystal::t_SmithTransform const transform = get_smith_transform( structure );
            Crystal::Structure::t_Atoms::iterator i_atom = structure.atoms.begin();
            Crystal::Structure::t_Atoms::iterator i_atom_end = structure.atoms.end();
            Eigen::Vector3i const &smith( boost::tuples::get<1>(transform) );
            for(; i_atom != i_atom_end; ++i_atom )
            {
              Eigen::Vector3i indices( Crystal::get_smith_index(transform, i_atom->pos) );
              const size_t index
              ( 
                indices[2] + smith(2) * ( indices[1] + smith(1) * indices[0] )
              );
              i_atom->type = line[index] == '0' ? -1.0: 1.0; 
            }
            _callable( structure );
            
            for(i_atom = structure.atoms.begin(); i_atom != i_atom_end; ++i_atom )
              i_atom->type = i_atom->type > 0e0 ? -1e0: 1e0;
            structure.name = '-' + structure.name;
            _callable( structure );
          }
          while( (not file.eof()) and file.good() );
        }
        catch( ... )
        {
          LADA_PYTHON_ERROR
          ( 
            PyExc_RuntimeError,
            ( "Error while enumerating gus-file " + _file + "\n" )
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
      boost::python::def
      ( 
        "enum_gusfile",
        &details::enumerate_gusfile,
        ( boost::python::arg("filename"), boost::python::arg("callable") ),
        "Reads a gus-file and calls a "
        "\"callable( const LaDa::Crystal::Structure )->types::t_real\" for each structure."
      );
    }

  }
} // namespace LaDa

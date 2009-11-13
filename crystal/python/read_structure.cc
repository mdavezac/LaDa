//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <sstream>
#include <iterator>
#include <set>
#include <algorithm>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/format.hpp>

#include <boost/python/str.hpp>
#include <boost/python/def.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/str.hpp>

#include <python/debug.hpp>
#include <physics/physics.h>

#include "../read_structure.h"
#include "../read_poscar.h"

namespace LaDa
{
  namespace Python
  {
    namespace bp = boost::python;
    namespace details
    {
      void print_poscar( Crystal::TStructure<std::string> const &_structure, 
                         boost::python::tuple const &_species,
                         const std::string &_path )
      {
        namespace bf = boost::filesystem;
        typedef Crystal::TStructure<std::string> t_Structure;
        // Extract specie types.
        std::vector<std::string> species;
        try 
        {
          for( size_t i(0); i < bp::len(_species); ++i )
          {
            std::string const type = bp::extract<std::string>( _species[i] );
            species.push_back( type );
            if( Physics::Atomic::Z( type ) == 0 )
            {
              LADA_PYTHON_ERROR
              ( 
                PyExc_RuntimeError, 
                "Atomic type " + type + " is unknown.\n" 
              );
              boost::python::throw_error_already_set();
              return;
            }
          }
        }
        catch(...)
        {
          LADA_PYTHON_ERROR( PyExc_RuntimeError, "Could not extract species.\n" );
          boost::python::throw_error_already_set();
          return;
        }
//       { // check all species are included
//         std::set<std::string> a;
//         foreach( t_Structure::t_Atom const &atom, _structure.atoms )
//           a.insert( atom.type );
//         std::set<std::string> b;
//         foreach( std::string const &sp, species ) b.insert( sp );
//         std::set<std::string> c;
//         std::set_difference( b.begin(), b.end(), a.begin(), a.end(),
//                              std::inserter( c, c.begin() ) );
//         if( not c.empty() )
//         {
//           std::string error = "Species in structure are not all included in specie list:";
//           foreach( std::string const d, c ) error += " " + d;
//           error += ".\n";
//           LADA_PYTHON_ERROR( PyExc_RuntimeError, error );
//           bp::throw_error_already_set();
//         }
//       }

        // Checks path
        std::string POSCAR("POSCAR");
        bf::path path( _path );
        if( _path.size() >= POSCAR.size() )
          if( _path.substr( _path.size() - POSCAR.size() ) == POSCAR )
            path = _path.substr(0, _path.size() - POSCAR.size() ); 
        if( path.string().size() != 0 ) 
        {
          if( not ( bf::exists(path) and bf::is_directory(path) ) )
          {
            LADA_PYTHON_ERROR( PyExc_RuntimeError, "Don't understand path name.\n" );
            boost::python::throw_error_already_set();
            return;
          }
        }

        // open poscar file
        path = path / "POSCAR";
        std::ofstream file;
        try { file.open( path.string().c_str(), std::ofstream::out|std::ofstream::trunc ); }
        catch(...)
        {
          LADA_PYTHON_ERROR( PyExc_RuntimeError, "Could not open file: " + path.string() + ".\n" );
          boost::python::throw_error_already_set();
          return;
        }

        // prints name, scale, and cell.
        file << boost::format("%s\n%.18f\n"
                              "  %25.18f  %25.18f  %25.18f  \n"
                              "  %25.18f  %25.18f  %25.18f  \n"
                              "  %25.18f  %25.18f  %25.18f  \n")
                  % _structure.name % _structure.scale
                  % _structure.cell(0,0) % _structure.cell(1,0) % _structure.cell(2, 0)
                  % _structure.cell(0,1) % _structure.cell(1,1) % _structure.cell(2, 1)
                  % _structure.cell(0,2) % _structure.cell(1,2) % _structure.cell(2, 2);

        // print nb per specie
        foreach( std::string const sp, species )
        {
          t_Structure::t_Atoms::const_iterator i_first( _structure.atoms.begin() );
          t_Structure::t_Atoms::const_iterator const i_end( _structure.atoms.end() );
          size_t N(0);
          for(; i_first != i_end; ++i_first )
            if( i_first->type == sp ) ++N;
          file << N << " ";
        }
        file << "\nDirect\n";
        // print cartesian coord.
        foreach( std::string const sp, species )
        {
          atat::rMatrix3d const inv( !_structure.cell );
          t_Structure::t_Atoms::const_iterator i_first( _structure.atoms.begin() );
          t_Structure::t_Atoms::const_iterator const i_end( _structure.atoms.end() );
          for( size_t N(0); i_first != i_end; ++i_first )
            if( i_first->type == sp )
            {
              atat::rVector3d const vec(inv * i_first->pos); 
              file << boost::format("  %25.18f %25.18f %25.18f\n") % vec(0) % vec(1) %vec(2);
            }
        }
        file.close();
      };

      template<class T_TYPE>
        boost::shared_ptr< Crystal::TStructure<T_TYPE> >
          read_poscar(bp::tuple const &_species, std::string const &_path, bool _c)
          {
            namespace bp = boost::python;
            typedef Crystal::TStructure<T_TYPE> t_Structure;
            try
            { 
              boost::shared_ptr<t_Structure> result( new t_Structure );
              std::vector<T_TYPE> species;
              for( size_t i(0); i < bp::len(_species); ++i )
              {
                T_TYPE const type = bp::extract<T_TYPE>( _species[i] );
                species.push_back( type );
              }
              Crystal::read_poscar( *result, _path, species, _c); 
              if( result->atoms.size() == 0 )
              {
                LADA_PYTHON_ERROR( PyExc_RuntimeError, "Could not read file " + _path + ".\n" );
                boost::python::throw_error_already_set();
                return boost::shared_ptr< Crystal::TStructure<T_TYPE> >(); 
              }
              return result;
            }
            catch( std::exception &_e )
            {
              LADA_PYTHON_ERROR( PyExc_RuntimeError, _e.what() );
              boost::python::throw_error_already_set();
              return boost::shared_ptr< Crystal::TStructure<T_TYPE> >(); 
            }
            catch( ... )
            {
              LADA_PYTHON_ERROR
              ( 
                PyExc_RuntimeError,
                ("Could not read MBCE-type structure from " + _path + ".\n")
              );
              boost::python::throw_error_already_set();
              return boost::shared_ptr< Crystal::TStructure<T_TYPE> >(); 
            }
          }

      template< class T_TYPE >
        boost::shared_ptr< Crystal::TStructure<T_TYPE> >
          read_poscar_one_specie(bp::str const &_specie, const std::string &_path, bool _c)
          {
            namespace bp = boost::python;
            typedef Crystal::TStructure<T_TYPE> t_Structure;
            try
            { 
              boost::shared_ptr<t_Structure> result( new t_Structure );
              T_TYPE const specie = bp::extract<T_TYPE>( _specie );
              std::vector<T_TYPE> const species(1, specie);
              Crystal::read_poscar( *result, _path, species, _c ); 
              if( result->atoms.size() == 0 )
              {
                LADA_PYTHON_ERROR( PyExc_RuntimeError, "Could not read file " + _path + ".\n" );
                boost::python::throw_error_already_set();
                return boost::shared_ptr< Crystal::TStructure<T_TYPE> >(); 
              }
              return result;
            }
            catch( std::exception &_e )
            {
              LADA_PYTHON_ERROR( PyExc_RuntimeError, _e.what() );
              boost::python::throw_error_already_set();
              return boost::shared_ptr< Crystal::TStructure<T_TYPE> >(); 
            }
            catch( ... )
            {
              LADA_PYTHON_ERROR
              ( 
                PyExc_RuntimeError,
                ("Could not read MBCE-type structure from " + _path + ".\n")
              );
              boost::python::throw_error_already_set();
              return boost::shared_ptr< Crystal::TStructure<T_TYPE> >(); 
            }
          }

      boost::shared_ptr<Crystal::Structure> read_structure(const std::string &_path, bool _c)
      {
        try
        { 
          boost::shared_ptr<Crystal::Structure> result( new Crystal::Structure );
          Crystal::read_structure(*result, _path, _c); 
          return result;
        }
        catch( std::exception &_e )
        {
          LADA_PYTHON_ERROR( PyExc_RuntimeError, _e.what() );
          boost::python::throw_error_already_set();
          return boost::shared_ptr<Crystal::Structure>(); 
        }
        catch( ... )
        {
          LADA_PYTHON_ERROR
          ( 
            PyExc_RuntimeError,
            ("Could not read MBCE-type structure from " + _path + ".\n")
          );
          boost::python::throw_error_already_set();
          return boost::shared_ptr<Crystal::Structure>(); 
        }
      }

      boost::shared_ptr<Crystal::Structure> read_pifile_structure( std::string const& _str )
      {
        try
        {
          std::istringstream sstring( _str );
          boost::shared_ptr<Crystal::Structure> structure( new Crystal::Structure );
          Crystal::read_pifile_structure( sstring, *structure );
          std::cout << _str << "\n";
          std::cout << *structure << "\n\n";
          return structure;
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
        return boost::shared_ptr<Crystal::Structure>();
      }
    }

    void expose_read_structure()
    {
      boost::python::def
      (
        "print_poscar", 
        &details::print_poscar,
        ( bp::arg("structure"), bp::arg("species"), bp::arg("path") = "./"),
        "Prints out a poscar to file." 
      );
      boost::python::def
      (
        "read_poscar", 
        &details::read_poscar_one_specie<std::string>,
        (
          bp::arg("types"),
          bp::arg("path"),
          bp::arg("check_lattice") = false
        ),
        "Tries to read a VASP POSCAR. Needs a tuple of species on input (first argument).\n"
        "Returns a structure on success" 
      );
      boost::python::def
      (
        "read_poscar", 
        &details::read_poscar<std::string>,
        (
          bp::arg("path"),
          bp::arg("types"),
          bp::arg("check_lattice") = false
        ),
        "Tries to read a VASP POSCAR. This is the case of a single atomic specie."
        "Returns a structure on success" 
      );
      boost::python::def
      (
        "read_structure", 
        &details::read_structure,
        (
          bp::arg("path"),
          bp::arg("check_lattice") = false
        ),
        "Tries to read a file as a MBCE structure file. Returns a structure on success" 
      );
      boost::python::def
      ( 
        "read_pifile_structure",
        &details::read_pifile_structure,
        bp::arg("string"),
        "Reads a structure from pi-file type input.\n" 
      );

    }

  }
} // namespace LaDa

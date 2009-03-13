//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>

#include <boost/python/self.hpp> 
#include <boost/python/class.hpp> 
#include <boost/python/def.hpp> 
#include <boost/python/list.hpp> 
#include <boost/python/str.hpp> 
#include <boost/python/extract.hpp> 
#include <boost/python/return_arg.hpp>
#include <boost/python/manage_new_object.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp> 
#include <boost/python/suite/indexing/vector_indexing_suite.hpp> 

#include <opt/types.h>
#include <python/std_map.hpp>
#include "../clj.h"

namespace LaDa
{
  namespace Python
  {
    template<class T_TYPE>
      void unfold_structure( const Crystal :: TStructure<T_TYPE>& _structure, std::vector<types::t_real>& _list )
      {
        _list.clear();
        typedef Crystal::TStructure<T_TYPE> t_Str;
        for( size_t i(0); i < 3; ++i )
          for( size_t j(0); j < 3; ++j )
            if( (_structure.freeze & Crystal::cell_freeze_tag(i,j)) == 0 )
              _list.push_back( _structure.cell(i,j) );
        typename Crystal::TStructure<T_TYPE>::t_Atoms::const_iterator i_atom = _structure.atoms.begin();
        typename Crystal::TStructure<T_TYPE>::t_Atoms::const_iterator i_atom_end = _structure.atoms.end();
        for(; i_atom != i_atom_end; ++i_atom )
        {
          if( (i_atom->freeze & t_Str::t_Atom::FREEZE_X) == 0 )
            _list.push_back( i_atom->pos[0] );
          if( (i_atom->freeze & t_Str::t_Atom::FREEZE_Y) == 0 )
            _list.push_back( i_atom->pos[1] );
          if( (i_atom->freeze & t_Str::t_Atom::FREEZE_Z) == 0 )
            _list.push_back( i_atom->pos[2] );
        }
      }
    template<class T_TYPE>
      void fold_structure( const std::vector<types::t_real> &_list,
                           Crystal :: TStructure<T_TYPE>& _structure )
      {
        namespace bp = boost::python;
        typedef Crystal::TStructure<T_TYPE> t_Str;
        std::vector<types::t_real> :: const_iterator i_var = _list.begin();
        for( size_t i(0); i < 3; ++i )
          for( size_t j(0); j < 3; ++j )
            if( (_structure.freeze & Crystal::cell_freeze_tag(i,j)) == 0 )
              { _structure.cell(i,j) = *i_var; ++i_var; }
        typename Crystal::TStructure<T_TYPE>::t_Atoms::iterator i_atom = _structure.atoms.begin();
        typename Crystal::TStructure<T_TYPE>::t_Atoms::iterator i_atom_end = _structure.atoms.end();
        for(; i_atom != i_atom_end; ++i_atom )
        {
          if( (i_atom->freeze & t_Str::t_Atom::FREEZE_X) == 0 )
            { i_atom->pos[0] = *i_var; ++i_var; }
          if( (i_atom->freeze & t_Str::t_Atom::FREEZE_Y) == 0 )
            { i_atom->pos[1] = *i_var; ++i_var; }
          if( (i_atom->freeze & t_Str::t_Atom::FREEZE_Z) == 0 )
            { i_atom->pos[2] = *i_var; ++i_var; }
        }
      }

    void read_fortran_input( Models::Clj &_clj, boost::python::list &_atoms, 
                             const boost::python::str& _path )
    {
      const std::string str = boost::python::extract<std::string>( _path );
      const boost::filesystem::path path( str );
      std::vector<std::string> atoms;
      Models::read_fortran_input( _clj, atoms, path );
      foreach( const std::string &atom, atoms )
        _atoms.append(atom);
    }

    template<class T_TYPE>
      void clear_structure( Crystal :: TStructure<T_TYPE>& _structure )
      {
        _structure.cell.zero();
        typename Crystal::TStructure<T_TYPE>::t_Atoms::iterator i_atom = _structure.atoms.begin();
        typename Crystal::TStructure<T_TYPE>::t_Atoms::iterator i_atom_end = _structure.atoms.end();
        for(; i_atom != i_atom_end; ++i_atom )
          for( size_t j(0); j < 3; ++j )
            i_atom->pos[j] = 0e0;
      }


    void expose_clj()
    {
      namespace bp = boost :: python;
      typedef Models::Clj::t_Bonds t_Bonds; 
      typedef Models::Clj::t_Bonds::value_type::second_type t_Bond; 
      typedef Models::Clj::t_Charges t_Charges; 

      expose_map< t_Charges >( "Test", "Dictionary of charges" );
      expose_map< t_Charges >( "Charges", "Dictionary of charges" );
//     bp::class_< t_Charges >( "Charges", "Dictionary of charges" )
//       .def( bp :: map_indexing_suite< t_Charges >() );
   
      bp::class_< t_Bond >( "LJBond", "Holds bond parameters for Lennard-Jones." )
        .def( bp::init<>() )
        .def( bp::init<t_Bond>() )
        .def_readwrite( "hardsphere", &t_Bond::hard_sphere, "+r^12 factor." )
        .def_readwrite( "vandderwalls", &t_Bond::van_der_walls, "-r^6 factor." );

      bp::class_< t_Bonds >( "LJBonds", "Dictionary of bonds" )
        .def( bp :: map_indexing_suite< t_Bonds >() );

      bp::class_< Models::Clj >( "Clj", "Coulomb + LennardJones functional.\n" )
        .def_readwrite( "charges", &Models::Clj::charges, "Dictionnary of charges." )
        .def_readwrite( "bonds", &Models::Clj::bonds, "Dictionnary of charges." )
        .def("__call__", &Models::Clj::operator() )
        .def("gradient", &Models::Clj::gradient );
 
      bp::def
      ( 
        "unfold_structure",
        &unfold_structure< std::string >,
        "Unfolds a structure into a c++ vector."
      );
      bp::def
      ( 
        "fold_structure",
        &fold_structure< std::string >,
        "Folds a c++ vector into a structure."
      );

      bp::def
      (
        "read_epinput",
        &read_fortran_input,
        (
          bp::arg("functional"), 
          bp::arg("species"),
          bp::arg("filename")
        ),
        "Reads input from fortran model.\n"
        "species is a list of atomic symbols which can be used to read a POSCAR."
      );
    }
  } // namespace Python
} // namespace LaDa

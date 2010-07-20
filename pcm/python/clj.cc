#include "LaDaConfig.h"

#include <vector>

#include <boost/python/self.hpp> 
#include <boost/python/class.hpp> 
#include <boost/python/def.hpp> 
#include <boost/python/list.hpp> 
#include <boost/python/str.hpp> 
#include <boost/python/extract.hpp> 
#include <boost/python/return_arg.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/manage_new_object.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp> 
#include <boost/python/suite/indexing/vector_indexing_suite.hpp> 
#include <boost/tuple/tuple.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/map.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <opt/types.h>
#include <python/misc.hpp>
#include "../clj.h"

namespace LaDa
{
  namespace python
  {
    namespace bp = boost::python;
    template<class T_TYPE>
      void unfold_structure( const Crystal :: TStructure<T_TYPE>& _structure,
                             std::vector<types::t_real>& _list )
      {
        _list.clear();
        typedef Crystal::TStructure<T_TYPE> t_Str;
        for( size_t i(0); i < 3; ++i )
          for( size_t j(0); j < 3; ++j )
            if( (_structure.freeze & Crystal::cell_freeze_tag(i,j)) == 0 )
              _list.push_back( _structure.cell(i,j) );
        typedef typename Crystal::TStructure<T_TYPE>::t_Atoms::const_iterator t_cit;
        t_cit i_atom = _structure.atoms.begin();
        t_cit i_atom_end = _structure.atoms.end();
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

    Models::Clj::t_Bonds::value_type::second_type* bond_constructor( types::t_real _a,
                                                                     types::t_real _b )
    {
      typedef Models::Clj::t_Bonds::value_type::second_type t_bond;
      t_bond *result( new t_bond );
      result->hard_sphere = _a;
      result->van_der_walls = _b;
      return result;
    }
    Models::Clj::t_Bonds::value_type::second_type* bond_tconstructor( bp::tuple const &_t )
    {
      if( bp::len( _t ) != 2 )
      {
        PyErr_SetString(PyExc_ValueError, "Expected a 3-tuple.");
        return NULL;
      }
      types::t_real const a = bp::extract<types::t_real>(_t[0]);
      types::t_real const b = bp::extract<types::t_real>(_t[1]);
      return bond_constructor(a, b);
    }

    void set_mesh( Models::Clj& _clj, const boost::python::tuple& _tuple )
    {
      namespace bp = boost::python;
      if( bp::len( _tuple ) != 3 )
      {
        PyErr_SetString(PyExc_ValueError, "Expected a 3-tuple.");
        return;
      }
      const types::t_int a = bp::extract<types::t_int>( _tuple[0] );
      const types::t_int b = bp::extract<types::t_int>( _tuple[1] );
      const types::t_int c = bp::extract<types::t_int>( _tuple[2] );
      _clj.set_mesh( boost::tuples::make_tuple( a, b, c ) );
    }
    boost::python::tuple get_mesh( Models::Clj& _clj )
    {
      namespace bt = boost::tuples;
      const boost::tuple<types::t_int, types::t_int, types::t_int> t( _clj.get_mesh() );
      return boost::python::make_tuple( bt::get<0>(t), bt::get<1>(t), bt::get<2>(t) ); 
    }

    template<int D> boost::shared_ptr<Models::Clj::t_Arg>
      __call__(Models::Clj const & _clj, Models::Clj::t_Arg const &_in)
      {
        boost::shared_ptr<Models::Clj::t_Arg> forces(new Models::Clj::t_Arg);
        forces->cell = math::rMatrix3d::Zero();
        forces->atoms.resize(_in.atoms.size());
        switch(D)
        {
          case 0: forces->energy = _clj.lennard_jones(_in, *forces); break;
          case 1: forces->energy = _clj.ewald(_in, *forces); break;
          case 2: forces->energy = _clj(_in, *forces); break;
        }
        return forces;
      }
    
    void set_sequence_item(Models::Clj::t_Bonds &_bonds, 
                           Models::Clj::t_Bonds::value_type::first_type const &_index,
                           bp::object const &_bond)
    {
      if(bp::len(_bond) != 2)
      {
        PyErr_SetString(PyExc_ValueError, "Incorrect value for bond type.");
        return;
      }
      types::t_real a, b;
      try
      { 
        a = bp::extract<types::t_real>(_bond[0]); 
        b = bp::extract<types::t_real>(_bond[1]); 
      }
      catch(...) 
      {
        PyErr_SetString(PyExc_ValueError, "Incorrect value for bond type.");
        return;
      }
      _bonds[_index] = Models::Clj::t_Bonds::value_type::second_type(a,b);
    }

    void expose_clj()
    {
      namespace bp = boost :: python;
      typedef Models::Clj::t_Bonds t_Bonds; 
      typedef Models::Clj::t_Bonds::value_type::second_type t_Bond; 
      typedef Models::Clj::t_Charges t_Charges; 

      bp::class_< t_Charges >( "Charges", "Dictionary of charges" )
        .def( bp :: map_indexing_suite< t_Charges >() )
        .def_pickle( Python::pickle< t_Charges >() );
   
      bp::class_< t_Bond >( "LJBond", "Holds bond parameters for Lennard-Jones." )
        .def( bp::init<>() )
        .def( bp::init<t_Bond>() )
        .def( "__init__", bp::make_constructor( &bond_constructor ),
              "First argument is hard-sphere parameter, and second vand-der-walls parameter." )
        .def( "__init__", bp::make_constructor( &bond_tconstructor ),
              "First argument is hard-sphere parameter, and second vand-der-walls parameter." )
        .def_readwrite( "hardsphere", &t_Bond::hard_sphere, "+r^12 factor." )
        .def_readwrite( "vanderwalls", &t_Bond::van_der_walls, "-r^6 factor." )
        .def_pickle( Python::pickle< t_Bond >() );

      bp::class_< t_Bonds >( "LJBonds", "Dictionary of bonds" )
        .def( bp :: map_indexing_suite< t_Bonds >() )
        .def("__setitem__", &set_sequence_item,
             " Sets lennard-jones parameter for a given bond.\n\n"
             " :Parameters:\n"
             "   - `index` of the bond, as obtained from `pcm.bond_name`.\n"
             "   - `value` a two tuple containing the facto of the r^12 term and the r^6 term.\n"
            )
        .def_pickle( Python::pickle< t_Bonds >() );

      bp::class_< Models::Clj >( "Clj", "Coulomb + LennardJones functional.\n" )
        .def_readwrite( "charges", &Models::Clj::charges, "Dictionnary of charges." )
        .def_readwrite( "bonds", &Models::Clj::bonds, "Dictionnary of charges." )
        .def( "lennard_jones", &__call__<0>, bp::arg("structure"),
              "Returns Lennard-Jones forces and energy.\n\n"
              ":Parameter:\n"
              " - `structure` Crystal structure for wich to compute energy.\n"
              ":return: a lada.crystal.Structure object containing forces and energy.\n" )
        .def( "ewald", &__call__<1>, bp::arg("structure"),
              "Returns Ewald forces and energy.\n\n"
              ":Parameter:\n"
              " - `structure` Crystal structure for wich to compute energy.\n"
              ":return: a lada.crystal.Structure object containing forces and energy.\n" )
        .def( "__call__", &__call__<2>, bp::arg("structure"),
              "Returns Ewald+Lennard-Jones forces and energy.\n\n"
              ":Parameter:\n"
              " - `structure` Crystal structure for wich to compute energy.\n"
              ":return: a lada.crystal.Structure object containing forces and energy.\n" )
        .def("__str__", &Python::tostream<Models::Clj> )
        .add_property
        (
          "mesh", 
          &get_mesh, &set_mesh,
          "Size of the superstructure for lennard-jones real-space sum."
        )
        .add_property
        (
          "ewald_cutoff", 
          &Models::Clj::Ewald::get_rcutoff, &Models::Clj::Ewald::set_rcutoff,
          "Sets cutoff for real-space ewald sum."
        )
        .add_property
        (
          "lj_cutoff", 
          &Models::Clj::LennardJones::get_rcutoff, &Models::Clj::LennardJones::set_rcutoff,
          "Sets cutoff for real-space lennard-jhones sum."
        )
        .def_pickle( Python::pickle< Models::Clj >() );
 
      bp::def( "bond_name", &Models::LennardJones::bondname );
      bp::def
      ( 
        "_unfold_structure",
        &unfold_structure< std::string >,
        "Unfolds a structure into a c++ vector."
      );
      bp::def
      ( 
        "_fold_structure",
        &fold_structure< std::string >,
        "Folds a c++ vector into a structure."
      );
    }
  } // namespace python
} // namespace LaDa

#include "LaDaConfig.h"

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/errors.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <sstream>

#include <opt/types.h>
#include <opt/debug.h>
#include <python/numpy_types.hpp>

#include "../traits.h"



#include "atom.hpp"


namespace LaDa
{
  namespace python
  {
    namespace bp = boost::python;
    typedef crystal::traits< std::set<std::string> >::t_Atom t_Atom;
    template<class T> 
      bp::object get_pos(typename crystal::traits<T>::t_Atom &_in)
        { return array_from_ptr<math::rVector3d::Scalar>(_in.pos.data(), 3); }
    template<class T> 
      void set_pos( typename crystal::traits<T>::t_Atom &_in,
                    boost::python::object const &_pos )
      {
        if(bp::len(_pos) != 3)
          PyException<error::input>::throw_error("Input object is not a sequence of 3.");
        _in.pos[0] = boost::python::exctract<math::rVector3d::Scalar>(_pos[0]);
        _in.pos[1] = boost::python::exctract<math::rVector3d::Scalar>(_pos[1]);
        _in.pos[2] = boost::python::exctract<math::rVector3d::Scalar>(_pos[2]);
      }

    template< class T_TYPE >
      void expose_typed_atom( const std::string &_name,
                              const std::string &_ds,
                              const std::string &_typeds )
      {
        namespace bp = boost::python;
        typedef typename crystal::traits< T_TYPE >::t_Atom t_Atom;
        bp::class_< t_Atom >
        ( 
          _name.c_str(), 
          ( 
              _ds 
            +  "\nThis object can be constructed from:\n\n"
               "- with no argument\n"
               "- another `" + _name + "` object (deepcopy)\n"
               "- a numpy 1d-vector of length 3, in which ``self.pos``" 
                 " is the only variable to be set.\n"
               "- same as above + a type value.\n"
               "- same as above + a site index.\n"
          ).c_str()
        ).def(bp::init<t_Atom const &>())
         .add_property
         (
           "pos", 
           bp::make_function(&get_pos<t_Atom>, bp::with_custodian_and_ward_postcall<1, 0>()),
           &set_pos<t_Atom>,
           "A 1-dimensional numpy array of length 3 containing atomic position in cartesian units."
         )
         .def_readwrite( "site",   &t_Atom::site,
                         "index of the \"site\" as referenced by a LaDa.Lattice object." )
         .def_readwrite( TypeAttr<T_TYPE>::name.c_str(),   &t_Atom::type, _typeds.c_str() )
         .def_readwrite( "freeze", &t_Atom::freeze )
         .def_pickle( python::pickle< t_Atom >() )
         .def( "__str__",  &print<t_Atom> );
      }

    void expose_atom()
    {
      namespace bp = boost::python;

      bp::enum_<AtomFreezeMixin::Mixin::frozen::type>( "FreezeAtom", "Tags to freeze atomic coordinates." )
        .value(       "none", AtomFreezeMixin::frozen::FREEZE_NONE )
        .value(          "x", AtomFreezeMixin::frozen::FREEZE_X )
        .value(          "y", AtomFreezeMixin::frozen::FREEZE_Y )
        .value(          "z", AtomFreezeMixin::frozen::FREEZE_Z )
        .value(       "type", AtomFreezeMixin::frozen::FREEZE_T )
        .value( "cartesians", AtomFreezeMixin::frozen::FREEZE_CARTESIANS )
        .value(        "all", AtomFreezeMixin::frozen::FREEZE_ALL )
        .export_values();

      expose_typed_atom<t_Atom::t_Type>
      (
        "Atom", 
        "Atom for which the type is specified as a set of strings.\n\n"
        "A set is a collection of unique strings.",
        "set of strings representing the atomic specie(s) at this site.\n\n"
        "This is a set in the sense of a collection of unique strings."
      );
    }
  }
} // namespace LaDa

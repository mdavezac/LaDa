#include "LaDaConfig.h"

#include <set>

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/errors.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/stl_iterator.hpp>

#include <opt/types.h>
#include <opt/debug.h>
#include <python/numpy_types.h>
#include <python/misc.hpp>
#include <python/raw_constructor.hpp>

#include "../atom.h"
#include "../traits.h"
#include "../is_container.h"



#include "atom.hpp"


namespace LaDa
{
  namespace python
  {
    namespace bp = boost::python;
    template<class T> struct atom
    {
      typedef typename crystal::traits::StructureData<T>::t_Atom type;
    };
    template<class T> 
      bp::object get_pos(typename atom<T>::type &_in)
        { return math::numpy::array_from_ptr<math::rVector3d::Scalar>(_in.pos.data(), 3); }
    template<class T> 
      void set_pos( typename crystal::traits::StructureData<T>::t_Atom &_in,
                    boost::python::object const &_pos )
      {
        if(bp::len(_pos) != 3)
          PyException<error::TypeError>::throw_error("Input object is not a sequence of 3.");
        _in.pos[0] = boost::python::extract<math::rVector3d::Scalar>(_pos[0]);
        _in.pos[1] = boost::python::extract<math::rVector3d::Scalar>(_pos[1]);
        _in.pos[2] = boost::python::extract<math::rVector3d::Scalar>(_pos[2]);
      }

    template<class T> typename boost::enable_if< crystal::details::is_scalar<T>, void >::type
      extract_type0(bp::object _in, T &_type)
        { _type = bp::extract<T>(_in); }
    template<class T> typename boost::enable_if< crystal::details::is_scalar<T> , void>::type
      extract_type1(bp::object _in, T &_type)
      {
        if(bp::len(_in) > 1) 
          PyException<error::TypeError>::throw_error("Atomic type needs be a scalar.");
        _type = bp::extract<T>(_in); 
      }
    template<class T> typename boost::enable_if< crystal::details::is_set<T>, void >::type
      extract_type0(bp::object _in, T &_type)
      {
        bp::stl_input_iterator<typename T::value_type> i_first(_in);
        bp::stl_input_iterator<typename T::value_type> const i_end;
        for(; i_first != i_end; ++i_first) _type.insert(*i_first);
      }
    template<class T> typename boost::enable_if< crystal::details::is_container<T>, void >::type
      extract_type0(bp::object _in, T &_type)
      {
        bp::stl_input_iterator<typename T::value_type> i_first(_in);
        bp::stl_input_iterator<typename T::value_type> const i_end;
        for(; i_first != i_end; ++i_first) _type.push_back(*i_first);
      }
    template<class T> typename boost::enable_if< crystal::details::is_iterable<T> , void>::type
      extract_type1(bp::object _in, T &_type)
        { extract_type0(_in, _type); }

    template<class T>
      std::string scalar_kind(typename atom<T>::type const &) { return "scalar"; }
    template<class T>
      std::string list_kind(typename atom<T>::type const &) { return "list"; }
    template<class T>
      std::string set_kind(typename atom<T>::type const &) { return "set"; }

    template<class T>
      std::auto_ptr<typename atom<T>::type> 
        constructor(bp::tuple _args, bp::dict _kwargs) 
        {
          typedef typename atom<T>::type t_Atom;
          std::auto_ptr<t_Atom> result(new t_Atom);
          if(bp::len(_args) == 0 and bp::len(_kwargs) == 0) return result;
          if(bp::len(_args) >= 3 and _kwargs.has_key("position"))
            PyException<error::TypeError>::throw_error
               ("Position provided both in tuple and keyword argument.");
          else if(bp::len(_args) >= 3 or _kwargs.has_key("position"))
          {
            bp::object pos;
            if(_kwargs.has_key("position")) pos = _kwargs.get("position");
            else pos = _args.slice(bp::slice_nil(), 3);
            if(bp::len(pos) != 3) 
              PyException<error::TypeError>::throw_error
                ("Position is not a list of three scalars.");
            result->pos[0] = boost::python::extract<math::rVector3d::Scalar>(pos[0]);
            result->pos[1] = boost::python::extract<math::rVector3d::Scalar>(pos[1]);
            result->pos[2] = boost::python::extract<math::rVector3d::Scalar>(pos[2]);
          }
          if(bp::len(_args) > 3 and _kwargs.has_key("type"))
            PyException<error::TypeError>::throw_error
               ("Types provided both in tuple and keyword argument.");
          else if(bp::len(_args) == 4 or _kwargs.has_key("type"))
          {
            bp::object type;
            if(_kwargs.has_key("type")) type = _kwargs.get("type");
            else type = _args[3];
            extract_type0(type, result->type);
          }
          else if(bp::len(_args) > 4)
            extract_type1(_args.slice(3, bp::slice_nil()), result->type);
          if(_kwargs.has_key("freeze"))
            result->freeze = bp::extract<int>(_kwargs.get("freeze"));
          if(_kwargs.has_key("site"))
            result->site = bp::extract<int>(_kwargs.get("site"));
          return result;
        }

    template< class T_TYPE >
      bp::class_<typename atom<T_TYPE>::type>
        expose_typed_atom( const std::string &_name,
                                const std::string &_ds,
                                const std::string &_typeds )
        {
          namespace bp = boost::python;
          typedef typename crystal::traits::StructureData< T_TYPE >::t_Atom t_Atom;
          return bp::class_< t_Atom >
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
           .def( "__init__", bp::raw_constructor(&constructor<T_TYPE>),
                 ( "Initialize an " + _name + ".\n\n"
                   ":Parameters:\n"
                   "  position : list\n    Atomic coordinates. "
                   "These quantities are accessible as a keyword or as the first "
                   "three arguments of the constructor.\n"
                   "  type\n    For AtomVec and AtomSet, this is a "
                   "list of string representing atomic species. For AtomStr, this "
                   "a single string. Can be accessed both as a keyword, or as the "
                   "4th argument (to nth argument, in case of AtomVec and "
                   "AtomSet) of the constructor.\n"
                   "  site : int\n    Site index. Can only be attained via a keyword only."
                   "  freeze : int\n    Site index. Can only be attained via a keyword only.").c_str())
           .add_property
           (
             "pos", 
             bp::make_function(&get_pos<T_TYPE>, bp::with_custodian_and_ward_postcall<1, 0>()),
             &set_pos<T_TYPE>,
             "A 1-dimensional numpy array of length 3 containing atomic position in cartesian units."
           )
           .def_readwrite( "site",   &t_Atom::site,
                           "index of the \"site\" as referenced by a LaDa.Lattice object." )
           .def_readwrite( "freeze", &t_Atom::freeze )
           .def_readwrite( "type", &t_Atom::type )
           .def_pickle( python::pickle< t_Atom >() )
           .def( "__str__",  &tostream<t_Atom> );
        }

    void expose_atom()
    {
      namespace bp = boost::python;
      import_array(); // needed for NumPy 

      bp::enum_<crystal::AtomFreezeMixin::frozen::type>
        ( "FreezeAtom", "Tags to freeze atomic coordinates." )
        .value(       "none", crystal::AtomFreezeMixin::frozen::NONE )
        .value(          "x", crystal::AtomFreezeMixin::frozen::X )
        .value(          "y", crystal::AtomFreezeMixin::frozen::Y )
        .value(          "z", crystal::AtomFreezeMixin::frozen::Z )
        .value(       "type", crystal::AtomFreezeMixin::frozen::T )
        .value( "cartesians", crystal::AtomFreezeMixin::frozen::CARTESIANS )
        .value(        "all", crystal::AtomFreezeMixin::frozen::ALL )
        .export_values();

      expose_typed_atom< std::set<std::string> >
      (
        "AtomSet", 
        "Atom for which the type is specified as a set of strings.\n\n"
        "A set is a collection of unique strings.",
        "set of strings representing the atomic specie(s) at this site.\n\n"
        "This is a set in the sense of a collection of unique strings."
      ).add_property("kind", &set_kind< std::set<std::string> >, "Occupation is a set of strings.");
      expose_typed_atom< std::vector<std::string> >
      (
        "AtomVec", 
        "Atom for which the type is specified as a list of strings.",
        "List of atomic species."
      ).add_property("kind", &list_kind< std::vector<std::string> >, "Occupation is a list of strings.");
      expose_typed_atom<std::string>
      (
        "AtomStr", 
        "Atom for which the type is specified as a strings.\n\n"
        "A set is a collection of unique strings.",
        "String representing the occupation."
      ).add_property("kind", &scalar_kind<std::string>, "Occupation is a string.");
    }
  }
} // namespace LaDa

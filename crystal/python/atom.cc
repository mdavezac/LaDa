#include "LaDaConfig.h"

#include <set>

#include <boost/get_pointer.hpp>
#include "../atom.h"
namespace boost
{
  namespace python { 
    template<class T> struct pointee< LaDa::crystal::Atom<T> > 
    {
      typedef LaDa::crystal::AtomData<T> type;
    };
namespace objects {
  inline LaDa::crystal::AtomData< std::set<std::string> >* get_pointer(LaDa::crystal::Atom< std::set<std::string> > const &_atom)
   { return _atom.get(); }
}}}


#include <boost/python/has_back_reference.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/errors.hpp>
#include <boost/python/reference_existing_object.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/stl_iterator.hpp>

#include <opt/types.h>
#include <opt/debug.h>
#include <python/numpy_types.h>
#include <python/misc.hpp>
#include <python/raw_constructor.hpp>


#include "atom.hpp"

namespace boost
{
  namespace python
  {

//     // So boost python knows how to call the constructors.
//     template<>
//       struct has_back_reference<typename LaDa::crystal::traits::StructureData<std::string>::t_Atom> : mpl::true_ {};
//     // So boost python knows how to call the constructors.
//     template<>
//       struct has_back_reference<
//         typename LaDa::crystal::traits::StructureData< std::vector<std::string> >::t_Atom> : mpl::true_ {};
//     // So boost python knows how to call the constructors.
//     template<>
//       struct has_back_reference<
//         typename LaDa::crystal::traits::StructureData< std::set<std::string> >::t_Atom> : mpl::true_ {};
  }
}

namespace LaDa
{
  namespace python
  {
    namespace bp = boost::python;
    template<class T> struct atom
    {
      typedef typename crystal::traits::StructureData<T>::t_Atom type;
      typedef typename type::element_type element;
    };
    template<class T> 
      bp::object get_pos(typename atom<T>::type &_in)
      {
        if(_in.self().ptr() == Py_None)
        {
           typename bp::reference_existing_object::apply<typename atom<T>::element *>::type converter;
           converter(_in.get());
//          self.attr("_set_self")();
        }
        return math::numpy::array_from_ptr<math::rVector3d::Scalar>(_in.pos.data(), 3, _in.self().ptr()); 
      }
    template<class T> 
      void set_pos( typename crystal::traits::StructureData<T>::t_Atom &_in,
                    boost::python::object _pos )
      {
        if(not is_position(_pos)) 
          PyException<error::TypeError>::throw_error("Input could not be converted to position.");
        extract_position(_pos, _in.pos);
      }

    template<class T>
      std::string scalar_kind(typename atom<T>::type const &) { return "scalar"; }
    template<class T>
      std::string list_kind(typename atom<T>::type const &) { return "list"; }
    template<class T>
      std::string set_kind(typename atom<T>::type const &) { return "set"; }

//   template<class T>
//     std::auto_ptr<typename atom<T>::type> constructor(bp::tuple _args, bp::dict _kwargs)
//     {
//       if(pb::len(_args)) == 0
//     }

//   template<class T>
//     std::auto_ptr<typename atom<T>::type> 
//       constructor(bp::tuple _args, bp::dict _kwargs) 
//       {
//         typedef typename atom<T>::type t_Atom;
//         std::auto_ptr<t_Atom> result(new t_Atom);
//         int current_index = 0;
//         if(bp::len(_args) == 0 and bp::len(_kwargs) == 0) return result;
//         else if(bp::len(_args) == 0) extract_atom(_kwargs, *result);
//         else if(is_position(_args[0]))
//         {
//           if(_kwargs.has_key("pos")) 
//             PyException<error::TypeError>("Constructor given position as argument and tuple.");
//           _kwargs["pos"] = _args[0];
//           current_index = 1;
//         }
//         else if(bp::len(_args) >= 3)
//         {
//           if(is_position(_args.slice(0, 3)))
//           {
//             if(_kwargs.has_key("pos")) 
//               PyException<error::TypeError>("Constructor given position as argument and tuple.");
//             _kwargs["pos"] = _args.slice(0, 3);
//             current_index = 3
//           }
//           else if(not _kwargs.has_key("pos"))
//             PyException<error::ValueError>::throw_error
//         }
//         else if(not _kwargs.has_key("pos"))
//           PyException<error::ValueError>::throw_error
//              ("Could not find position in arguments nor in keyword arguments.");
//         else if(bp::len(_args) == 4 or _kwargs.has_key("type"))
//         {
//           bp::object type;
//           if(_kwargs.has_key("type")) type = _kwargs.get("type");
//           else type = _args[3];
//           extract_specie(type, result->type);
//         }
//         else if(bp::len(_args) > 4)
//           details::extract_type1(_args.slice(3, bp::slice_nil()), result->type);
//         if(_kwargs.has_key("freeze"))
//           result->freeze = bp::extract<int>(_kwargs.get("freeze"));
//         if(_kwargs.has_key("site"))
//           result->site = bp::extract<int>(_kwargs.get("site"));
//         bp::stl_iterator<std::string> i_key(_kwargs.keys());
//         bp::stl_iterator<std::string> const i_key_end;
//         for(size_t i(0); i_key != i_key_end; ++i_key)
//           if(*i_key != "pos" and *i_key != "site" and *i_key != "freeze" and *i_key != "type")
//             result->__dict__[*i_key] = _kwargs.get(*i_key);
//
//         return result;
//       }

    template<class T>
      void set_self(bp::back_reference<typename atom<T>::type&> _ref)
      {
        if(_ref.get().self().ptr() == Py_None)
          _ref.get().set_self(_ref.source().ptr());
        else if(_ref.get().self().ptr() != _ref.source().ptr())
          PyException<error::InternalError>::throw_error("back_reference and self do not point to same object.");
      }

    template< class T_TYPE >
      bp::class_<typename atom<T_TYPE>::element, typename atom<T_TYPE>::type, boost::noncopyable>
        expose_typed_atom( const std::string &_name,
                                const std::string &_ds,
                                const std::string &_typeds )
        {
          namespace bp = boost::python;
          typedef typename crystal::traits::StructureData< T_TYPE >::t_Atom t_Atom;
          typedef typename t_Atom::element_type t_Element;
          return bp::class_<t_Element, t_Atom, boost::noncopyable>
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
          ); //.def(bp::init<t_Atom const &>())
          //.def( "__init__", bp::raw_constructor(&constructor<T_TYPE>),
          //      ( "Initialize an " + _name + ".\n\n"
          //        ":Parameters:\n"
          //        "  pos : list\n    Atomic coordinates. "
          //        "These quantities are accessible as a keyword or as the first "
          //        "three arguments of the constructor.\n"
          //        "  type\n    For AtomVec and AtomSet, this is a "
          //        "list of string representing atomic species. For AtomStr, this "
          //        "a single string. Can be accessed both as a keyword, or as the "
          //        "4th argument (to nth argument, in case of AtomVec and "
          //        "AtomSet) of the constructor.\n"
          //        "  site : int\n    Site index. Can only be attained via a keyword only."
          //        "  freeze : int\n    Site index. Can only be attained via a keyword only.").c_str())
//          .add_property
//            (
//              "pos", 
//              bp::make_function(&get_pos<T_TYPE>, bp::with_custodian_and_ward_postcall<1, 0>()),
//              &set_pos<T_TYPE>,
//              "A 1-dimensional numpy array of length 3 containing atomic position in cartesian units."
//            )
//          .def_readwrite( "site",   &t_Element::site,
//                          "index of the \"site\" as referenced by a LaDa.Lattice object." )
//          .def_readwrite( "freeze", &t_Element::freeze )
//          .def_readwrite( "type", &t_Element::type );
//          .def_pickle( python::pickle< t_Element >() )
//          .def("_set_self", &set_self<T_TYPE>)
//          .def( "__str__",  &tostream<t_Element> );
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
      ); //.add_property("kind", &set_kind< std::set<std::string> >, "Occupation is a set of strings.");
//     expose_typed_atom< std::vector<std::string> >
//     (
//       "AtomVec", 
//       "Atom for which the type is specified as a list of strings.",
//       "List of atomic species."
//     ).add_property("kind", &list_kind< std::vector<std::string> >, "Occupation is a list of strings.");
//     expose_typed_atom<std::string>
//     (
//       "AtomStr", 
//       "Atom for which the type is specified as a strings.\n\n"
//       "A set is a collection of unique strings.",
//       "String representing the occupation."
//     ).add_property("kind", &scalar_kind<std::string>, "Occupation is a string.");
//     bp::register_ptr_to_python< std::auto_ptr<atom<std::string>::type> >();
    }
  }
} // namespace LaDa

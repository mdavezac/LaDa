#ifndef LADA_CRYSTAL_PYTHON_ATOM_HPP
#define LADA_CRYSTAL_PYTHON_ATOM_HPP

#include "LaDaConfig.h"

#include <math/python/python.hpp>

#include "../traits.h"
#include "../is_container.h"

namespace LaDa
{
  namespace python
  {
    //! Declares interface to atom.
    void expose_atom();

//   //! checks whether this is a specie index.
//   template<class T> typename boost::enable_if< crystal::details::is_scalar<T>, bool > :: type
//     is_specie(bp::object const &_index)
//       { return bp::extract<T>(_index).check(); }
//   //! checks whether this is a specie index.
//   template<class T> typename boost::enable_if< crystal::details::is_iterable<T>, bool > :: type
//     is_specie(bp::object const &_index)
//     {
//       if(bp::extract<T&>(_index).check()) return true;
//       if(bp::extract<typename T::value_type&>(_index).check()) return true;
//       if(not PySequence_Check(_index.ptr()))
//         PyException<error::TypeError>::throw_error("Input is not a sequence.");
//       bp::object i_specie( bp::handle<>(PyObject_GetIter(_index.ptr())) );
//       if(i_specie.ptr() == Py_None) return true;
//       while(PyObject *specie = PyIter_Next(i_specie.ptr()))
//       {
//         if(not bp::extract<typename T::value_type>(specie).check())
//         {
//           Py_DECREF(specie);
//           return false;
//         }
//         Py_DECREF(specie);
//       }
//       return true;
//     };
//   //! Extracts scalar type to c++.
//   template<class T> typename boost::enable_if< crystal::details::is_scalar<T>, void >::type
//     extract_specie(bp::object _in, T &_type)
//       { _type = bp::extract<T>(_in); }
//   //! Extracts set type to c++.
//   template<class T> typename boost::enable_if< crystal::details::is_set<T>, void >::type
//     extract_specie(bp::object _in, T &_type)
//     {
//       if(bp::extract<T&>(_in).check())
//         _type = bp::extract<T&>(_in);
//       else if(bp::extract<typename T::value_type&>(_index).check())
//         _type.insert( bp::extract<typename T::value_type const &>(_in)() );
//       else
//       {
//         bp::stl_input_iterator<typename T::value_type> i_first(_in);
//         bp::stl_input_iterator<typename T::value_type> const i_end;
//         for(; i_first != i_end; ++i_first) _type.insert(*i_first);
//       }
//     }
//   //! Extracts vector type to c++.
//   template<class T> typename boost::enable_if< crystal::details::is_container<T>, void >::type
//     extract_specie(bp::object _in, T &_type)
//     {
//       if(bp::extract<T&>(_in).check())
//         _type = bp::extract<T&>(_in);
//       else if(bp::extract<typename T::value_type&>(_index).check())
//         _type.push_back( bp::extract<typename T::value_type&>(_in)() );
//       else
//       {
//         bp::stl_input_iterator<typename T::value_type> i_first(_in);
//         bp::stl_input_iterator<typename T::value_type> const i_end;
//         for(; i_first != i_end; ++i_first) _type.push_back(*i_first);
//       }
//     }
//   
//   //! checks whether this is convertible to an atom.
//   template<class T> 
//     bool is_atom(bp::object const &_index)
//     {
//       bp::object index(_index);
//       // Actual atom object.
//       if(bp::extract<typename crystal::TemplateStructure<T>::reference>(_index).check())
//         return true;
//       // Object as dictionary.
//       else if(PyDict_Check(index.ptr()))
//       {
//         bp::dict const dict(_index);
//         if(not dict.has_key("pos")) return false;
//         if(not is_position(dict["pos"])) return false;
//         if(not dict.has_key("type")) return false;
//         if(not is_specie<T>(dict["type"])) return false;
//         return true;
//       }
//       // Object as sequence.
//       else if(PySequence_Check(index.ptr()))
//       {
//         if(bp::len(_index) >= 4) return false;
//         if(not bp::extract<math::rMatrix3d::Scalar>(_index[0]).check()) return false;
//         if(not bp::extract<math::rMatrix3d::Scalar>(_index[1]).check()) return false;
//         if(not bp::extract<math::rMatrix3d::Scalar>(_index[2]).check()) return false;
//         if(is_specie<T>(_index[3])) return bp::len(_index) == 4;
//         if(crystal::details::is_scalar<T>::value) return false;
//         return is_specie<T>(_index.slice(3, bp::slice_nil));
//       }
//       // Object as object.
//       if(not PyObject_HasAttrString(index.ptr(), "pos")) return false;
//       if(not is_position(index.attr("pos"))) return false;
//       if(not PyObject_HasAttrString(index.ptr(), "type")) return false;
//       return is_specie<T>(index.attr("type"))
//     }
//
//   //! Creates an atom out of a dictionary, sequence, or object.
//   template<class T> 
//     void extract_atom(bp::object const _in, crystal::TemplateStructure<T>::reference _atom)
//     {
//       if(bp::extract<typename crystal::TemplateStructure<T>::reference>(_index).check())
//         _atom = bp::extract<typename crystal::TemplateStructure<T>::const_reference>(_index);
//       else if(PyDict_Check(_in.ptr()))
//       {
//         bp::dict const dict(bp::handle<>(PyDict_Copy(_in.ptr())));
//         extract_position(dict["pos"], _atom.pos);
//         PyDict_DelItemString(dict.ptr(), "pos");
//         extract_specie(dict["type"], _atom.type);
//         PyDict_DelItemString(dict.ptr(), "type");
//         if(dict.has_key("site"))
//         {
//           _atom.site = bp::extract<types::t_int>(dict["site"]);
//           PyDict_DelItemString(dict.ptr(), "site");
//         }
//         if(dict.has_key("freeze"))
//         {
//           _atom.site = bp::extract<types::t_int>(dict["freeze"]);
//           PyDict_DelItemString(dict.ptr(), "freeze");
//         }
//         bp::dict nmspace; nmspace["input"] = dict;
//         bp::dict copy_dict( bp::handle<> (
//           "from copy import deepcopy\n"
//           "result = deepcopy(input)",
//           Py_file_input,
//           nmspace.ptr(),
//           nmspace.ptr()
//           ) );
//         copy_dict = copy_dict["input"];
//         bp::stl_iterator<std::string> i_key(_kwargs.keys());
//         bp::stl_iterator<std::string> const i_key_end;
//         for(size_t i(0); i_key != i_key_end; ++i_key)
//           _atom.__dict__[*i_key] = copy_dict[*i_key];
//         return;
//       }
//       else if(PySequence_Check(_in.ptr()))
//       {
//         extract_position(_in.slice(0, 3));
//         if(crystal::details::is_scalar<T>::value or bp::len(_in) == 4)
//           extract_specie<T>(_in[3], _atom.type);
//         else extract_specie<T>(_in.slice(3, bp::slice_nil), _atom.type);
//         return;
//       };
//       extract_position(_in.attr("pos"), _atom.pos);
//       extract_specie(_in.attr("type"), _atom.type);
//       if(PyObject_HasAttrString(_in.ptr(), "site"))
//         _atom.site = bp::extract<types::t_int>(_in.attr("site"));
//       if(PyObject_HasAttrString(_in.ptr(), "freeze"))
//         _atom.freeze = bp::extract<types::t_int>(_in.attr("freeze"));
//     };
//
//   //! Implementation details of python stuff.
//   namespace details
//   {
//     template<class T> typename boost::enable_if< crystal::details::is_scalar<T> , void>::type
//       extract_type1(bp::object _in, T &_type)
//       {
//         if(bp::len(_in) > 1) 
//           PyException<error::TypeError>::throw_error("Atomic type needs be a scalar.");
//         _type = bp::extract<T>(_in); 
//       }
//     template<class T> typename boost::enable_if< crystal::details::is_iterable<T> , void>::type
//       extract_type1(bp::object _in, T &_type)
//         { extract_type0(_in, _type); }
//   }
  }
} 
#endif

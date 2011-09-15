#include "LaDaConfig.h"

#include <set>
#include <algorithm>
#include <iterator>
#include <memory>

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/manage_new_object.hpp>
#include <boost/python/stl_iterator.hpp>

#include <opt/types.h>
#include <opt/debug.h>
#include <python/numpy_types.h>
#include <python/exceptions.h>

#include "../traits.h"

#include "set.hpp"


namespace LaDa
{
  namespace python
  {
    typedef std::set<std::string> type;
    namespace bp = boost::python;
    bool contains(type const &_in, type::value_type _key)
      { return _in.count(_key) == size_t(1); }

    bool isdisjoint(type const &_a, bp::object _b)
    {
      if(not PyObject_HasAttrString(_b.ptr(), "__iter__"))
        python::PyException<error::ArgumentError>::throw_error("Object is not a sequence.");
      bp::stl_input_iterator<type::value_type> i_first(_b);
      bp::stl_input_iterator<type::value_type> const i_end;
      for(; i_first != i_end; ++i_first)
        if(_a.count(*i_first) == 1) return false;
      return true;
    }

    bool issubset(type const &_a, bp::object _b)
    {
      if(PyObject_HasAttrString(_b.ptr(), "__contains__") == 0)
        python::PyException<error::ArgumentError>::throw_error("No __contains__ protocol.");
      bp::object contains = _b.attr("__contains__");
      type::const_iterator i_first = _a.begin();
      type::const_iterator const i_end = _a.end();
      for(; i_first != i_end; ++i_first)
        if(not bp::extract<bool>(contains(*i_first))) return false;
      return true;
    }
    bool istruesubset(type const &_a, bp::object _b)
    {
      if(PyObject_HasAttrString(_b.ptr(), "__len__") == 0)
        python::PyException<error::ArgumentError>::throw_error("No __len__ protocol.");
      if( _a.size() >= bp::len(_b) ) return false;
      return issubset(_a, _b);
    }
    bool issuperset(type const &_a, bp::object _b)
    {
      if(not PyObject_HasAttrString(_b.ptr(), "__iter__"))
        python::PyException<error::ArgumentError>::throw_error("No __iter__ protocol.");
      bp::stl_input_iterator<type::value_type> i_first(_b);
      bp::stl_input_iterator<type::value_type> const i_end;
      for(; i_first != i_end; ++i_first)
        if(_a.count(*i_first) == 0) return false;
      return true;
    }
    bool istruesuperset(type const &_a, bp::object _b)
    {
      if(PyObject_HasAttrString(_b.ptr(), "__len__") == 0)
        python::PyException<error::ArgumentError>::throw_error("No __len__ protocol.");
      if( _a.size() <= bp::len(_b) ) return false;
      return issuperset(_a, _b);
    }

    void insert_(type &_inout, bp::object _b)
    {
      if(not PyObject_HasAttrString(_b.ptr(), "__iter__"))
        python::PyException<error::ArgumentError>::throw_error("No __iter__ protocol.");
      bp::stl_input_iterator<type::value_type> i_first(_b);
      bp::stl_input_iterator<type::value_type> const i_end;
      for(; i_first != i_end; ++i_first) _inout.insert(*i_first);
    }

    bp::object union_(type const &_a, bp::object _b)
    {
      type result = _a;
      insert_(result, _b);
      return bp::object(result);
    }
    void update(type &_a, bp::object _b)
      { insert_(_a, _b); }

    bp::object intersection(type const &_a, bp::object _b)
    {
      type b;
      insert_(b, _b);
      type result;
      std::insert_iterator<type> inserter(result, result.begin());
      std::set_intersection(_a.begin(), _a.end(), b.begin(), b.end(), inserter);
      return bp::object(result);
    }

    void isub(type &_a, bp::object _b)
    {
      if(not PyObject_HasAttrString(_b.ptr(), "__iter__"))
        python::PyException<error::ArgumentError>::throw_error("No __iter__ protocol.");
      bp::stl_input_iterator<type::value_type> i_first(_b);
      bp::stl_input_iterator<type::value_type> const i_end;
      for(; i_first != i_end; ++i_first)
      {
        type::iterator i_found = _a.find(*i_first);
        if(i_found != _a.end()) _a.erase(i_found);
      }
    }

    bp::object difference(type const &_a, bp::object _b)
    {
      type result = _a;
      isub(result, _b);
      return bp::object(result);
    }

    bp::object symmetric_difference(type const &_a, bp::object _b)
    {
      type b;
      insert_(b, _b);
      type result;
      std::insert_iterator<type> inserter(result, result.begin());
      std::set_symmetric_difference(_a.begin(), _a.end(), b.begin(), b.end(), inserter);
      return bp::object(result);
    }
    void ixor(type &_a, bp::object _b)
    {
      type b;
      insert_(b, _b);
      type result;
      std::insert_iterator<type> inserter(result, result.begin());
      std::set_symmetric_difference(_a.begin(), _a.end(), b.begin(), b.end(), inserter);
      _a = result;
    }
    bp::object iand(type &_a, bp::object _b)
    {
      type b;
      insert_(b, _b);
      type result;
      std::insert_iterator<type> inserter(result, result.begin());
      std::set_intersection(_a.begin(), _a.end(), b.begin(), b.end(), inserter);
      _a = result;
      return bp::object(result);
    }

    bp::object copy(type const &_a) { return bp::object(type(_a)); }
    void add(type &_a, type::value_type const &_key) { _a.insert(_key); }
    void remove(type &_a, type::value_type const &_key)
    {
      type::iterator i_found = _a.find(_key);
      if(i_found == _a.end())
        python::PyException<error::KeyError>::throw_error("Atomic specie not in set.");
      _a.erase(i_found);
    }
    type::value_type pop(type &_a, type::value_type const &_key)
    {
      type::iterator i_found = _a.find(_key);
      if(i_found == _a.end())
        python::PyException<error::KeyError>::throw_error("Atomic specie not in set.");
      return _key;
    }
    void discard(type &_a, type::value_type const &_key)
    {
      type::iterator i_found = _a.find(_key);
      if(i_found != _a.end()) _a.erase(i_found);
    }

    std::auto_ptr<type> from_list(bp::object _in)
    {
      std::auto_ptr<type> result(new type);
      if(_in.ptr() == Py_None) return result;
      insert_(*result, _in);
      return result;
    }

    std::string __str__(type const &_a)
    {
      if( _a.size() == 0 ) return "SpecieSet()";
      std::ostringstream sstr;
      type::const_iterator i_first = _a.begin();
      type::const_iterator const i_end = _a.end();
      sstr << "SpecieSet([\"" << *i_first << "\"";
      for(++i_first; i_first != i_end; ++i_first)
        sstr << ", \"" << *i_first << "\"";
      sstr << "])";
      return sstr.str();
    }
    bool equal(type const &_a, bp::object _b)
    {
      type result;
      insert_(result, _b);
      return result == _a; 
    }

    void expose_set()
    {
      bp::class_<type>("SpecieSet", "Set of unique strings.")
        .def("__init__", bp::make_constructor(&from_list))
        .def("__len__", &type::size)
        .def("isdisjoint", &isdisjoint)
        .def("issubset", &issubset)
        .def("__le__", &issubset)
        .def("__lt__", &istruesubset)
        .def("issuperset", &issuperset)
        .def("__ge__", &issuperset)
        .def("__gt__", &istruesuperset)
        .def("union", &union_)
        .def("__or__", &union_)
        .def("intersection", &intersection)
        .def("__and__", &intersection)
        .def("difference", &difference)
        .def("__sub__", &difference)
        .def("symmetric_difference", &symmetric_difference)
        .def("__xor__", &symmetric_difference)
        .def("update", &update)
        .def("__ior__", &update)
        .def("intersection_update", &iand)
        .def("__iand__", &iand)
        .def("difference_update", &isub)
        .def("__isub__", &isub)
        .def("symmetric_difference_update", &ixor)
        .def("__ixor__", &ixor)
        .def("copy", &copy)
        .def("add", &add)
        .def("remove", &remove)
        .def("discard", &discard)
        .def("pop", &pop)
        .def("clear", &type::clear)
        .def("__eq__", &equal)
        .def("__iter__", bp::iterator<type>())
        .def("__repr__", &__str__)
        .def("__str__", &__str__);
    };
  }
}

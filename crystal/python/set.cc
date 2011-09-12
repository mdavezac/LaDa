#include "LaDaConfig.h"

#include <sstream>
#include <algorithm>
#include <iterator>

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/errors.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

#include <opt/types.h>
#include <opt/debug.h>
#include <python/numpy_types.hpp>

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

    bool isdisjoint(type const &_a, type const &_b)
    {
      type::const_iterator i_first = _b.begin();
      type::const_iterator const i_end = _b.end();
      for(; i_first != i_end; ++i_first)
        if(_a.count(*i_first) == 1) return false;
      return true;
    }

    bool issubset(type const &_a, type const &_b)
    {
      if(_a.size() > _b.size()) return false;
      type::const_iterator i_first = _a.begin();
      type::const_iterator const i_end = _a.end();
      for(; i_first != i_end; ++i_first)
        if(_b.count(*i_first) == 0) return false;
      return true;
    }
    bool istruesubset(type const &_a, type const &_b)
    {
      if(_b.size() == _a.size()) return false;
      return issubset(_a, _b);
    }
    bool issuperset(type const &_a, type const &_b)
      { return issubset(_b, _a); }
    bool istruesuperset(type const &_a, type const &_b)
      { return istruesubset(_b, _a); }

    bp::object union_(type const &_a, type const &_b)
    {
      type result;
      std::insert_iterator<type> inserter(result, result.begin());
      std::set_union(_a.begin(), _a.end(), _b.begin(), _b.end(), inserter);
      return bp::object(result);
    }

    bp::object intersection(type const &_a, type const &_b)
    {
      type result;
      std::insert_iterator<type> inserter(result, result.begin());
      std::set_intersection(_a.begin(), _a.end(), _b.begin(), _b.end(), inserter);
      return bp::object(result);
    }

    bp::object difference(type const &_a, type const &_b)
    {
      type result;
      std::insert_iterator<type> inserter(result, result.begin());
      std::set_difference(_a.begin(), _a.end(), _b.begin(), _b.end(), inserter);
      return bp::object(result);
    }

    bp::object symmetric_difference(type const &_a, type const &_b)
    {
      type result;
      std::insert_iterator<type> inserter(result, result.begin());
      std::set_symmetric_difference(_a.begin(), _a.end(), _b.begin(), _b.end(), inserter);
      return bp::object(result);
    }

    bp::object copy(type const &_a) { return bp::object(type(_a)); }

    void symmetric_difference(type const &_a, type const &_b)

    expose_set()
    {
      bp::class_<type>("SpecieSet", bp::no_init(), "Set of unique strings.")
        .def(bp::init<type const & >())
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
        .def("symmetric_difference", &symdiff)
        .def("__xor__", &symdiff)
        .def("copy", &copy)
    };
  }
}

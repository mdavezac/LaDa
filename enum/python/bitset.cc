//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <sstream>
#include <complex>

#include <boost/exception/diagnostic_information.hpp>

#include <boost/python/class.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/def.hpp>
#include <boost/python/register_ptr_to_python.hpp>

#include <python/std_vector.hpp>
#include <crystal/lattice.h>

#include "bitset.hpp"
#include "../numeric_type.h"



namespace LaDa
{
  
  namespace Python
  {
    namespace bp = boost::python;
    bool get_item( enumeration::Database const &_d, enumeration::t_uint _x ) { return _d[_x]; }
    void set_item( enumeration::Database &_d, enumeration::t_uint _x, bool _b ) { _d[_x] = _b; }
    void expose_bitset()
    {
      bp::def("get_index", &enumeration::get_index);
      bp::def("count_flavors", &enumeration::count_flavors);
      bp::def("create_flavorbase", &enumeration::create_flavor_base);
      expose_vector<size_t>("FlavorBase", "A basis k^m");
      bp::register_ptr_to_python< boost::shared_ptr< std::vector<size_t> > >();
      

      bp::scope scope = bp::class_<enumeration::Database>
      (
        "Database", 
        "A large bitset"
      ).def( bp::init<enumeration::Database const&>() )
       .def("__getitem__", &get_item)
       .def("__setitem__", &set_item)
       .def("__len__", &enumeration::Database::size);

      bp::def("__init__", &enumeration::create_database);
    }
  }
} // namespace LaDa

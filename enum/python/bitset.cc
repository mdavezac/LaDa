//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <sstream>
#include <complex>

#include <pyublas/numpy.hpp>

#include <boost/exception/diagnostic_information.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_signed.hpp>
#include <boost/type_traits/make_signed.hpp>

#include <boost/python/class.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/def.hpp>
#include <boost/python/object.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <boost/python/make_constructor.hpp>

#include <python/std_vector.hpp>
#include <crystal/lattice.h>

#include "bitset.hpp"
#include "../numeric_type.h"
#include "../exceptions.h"



namespace LaDa
{
  
  namespace Python
  {
    namespace bp = boost::python;
    bool get_item( enumeration::Database const &_d, enumeration::t_uint _x ) { return _d[_x]; }
    void set_item( enumeration::Database &_d, enumeration::t_uint _x, bool _b ) { _d[_x] = _b; }

    pyublas::numpy_vector<size_t> integer_to_vector( enumeration::t_uint _x,
                                                     enumeration::FlavorBase const& _fl ) 
    {
      pyublas::numpy_vector<size_t> result;
      enumeration::integer_to_vector(_x, _fl, result);
      return result;
    }

    std::string integer_to_bitstring(enumeration::t_uint _x, enumeration::FlavorBase const &_fl )
    {
      try
      {
        return enumeration::integer_to_bitstring(_x, _fl);
      }
      catch(enumeration::integer_too_large &_e)
      {
        std::ostringstream sstr;
        sstr << "Integer is too large according to FlavorBase.\n"
             << "Maximum size is " << _fl.back() * _fl[1] << "\n"
             << "Argument is " << _x << "\n"
             << boost::diagnostic_information(_e) << "\n";
        PyErr_SetString(PyExc_OverflowError, sstr.str().c_str());
        bp::throw_error_already_set();
        return "";
      }
      catch(boost::exception &_e)
      {
        std::ostringstream sstr;
        sstr << "Error from boost library.\n"
             << boost::diagnostic_information(_e);
        PyErr_SetString(PyExc_RuntimeError, sstr.str().c_str());
        bp::throw_error_already_set();
        return "";
      }
      catch(...)
      {
        PyErr_SetString(PyExc_ValueError, "Unknown error in integer_to_bitstring.\n");
        bp::throw_error_already_set();
        return "";
      }
    }

    boost::shared_ptr<enumeration::Database> create(size_t const &_a, size_t const &_b)
    {
      if( _a == 0 or _b == 0 )
      {
        PyErr_SetString(PyExc_ValueError, "Invalid argument: value of 0.\n" );
        bp::throw_error_already_set();
        return boost::shared_ptr<enumeration::Database>();
      }
      try { return enumeration::create_database(_a, _b); }
      catch(enumeration::supercell_too_large &_e)
      {
        std::ostringstream sstr;
        sstr << "Supercell configuration-space too large.\n"
             << "Maximum size is 2^" << sizeof(enumeration::t_uint) << "\n"
             << boost::diagnostic_information(_e);
        PyErr_SetString(PyExc_RuntimeError, sstr.str().c_str());
        bp::throw_error_already_set();
        return boost::shared_ptr<enumeration::Database>();
      }
      catch(boost::exception &_e)
      {
        std::ostringstream sstr;
        sstr << "Error from boost library.\n"
             << boost::diagnostic_information(_e);
        PyErr_SetString(PyExc_RuntimeError, sstr.str().c_str());
        bp::throw_error_already_set();
        return boost::shared_ptr<enumeration::Database>();
      }
      catch(...)
      {
        PyErr_SetString(PyExc_ValueError, "Unknown error while creating database.\n");
        bp::throw_error_already_set();
        return boost::shared_ptr<enumeration::Database>();
      }
    }

    template<class T> typename boost::enable_if<boost::is_signed<T>, bool>::type
      negative( T const &_a ) { return _a < 0; } 
    template<class T> typename boost::disable_if<boost::is_signed<T>, bool>::type
      negative( T const &_a ) { return false; }
    template<class T> typename boost::enable_if<boost::is_signed<T>, enumeration::t_uint>::type
      abs( T const &_a ) { return std::abs(_a); }
    template<class T> typename boost::disable_if<boost::is_signed<T>, enumeration::t_uint>::type
      abs( T const &_a ) { return _a; }

    template< class T1, class T2>
      boost::shared_ptr<enumeration::Database> create1(T1 const &_a, T2 const &_b)
      {
        if( negative(_a) or negative(_b) )
        {
          PyErr_SetString(PyExc_ValueError, "Invalid argument: negative value.\n" );
          bp::throw_error_already_set();
          return boost::shared_ptr<enumeration::Database>();
        }
        return create(_a, _b);
      }


    template<class T_STR> 
      inline void as_structure( T_STR &_out, enumeration::t_uint _x,
                                enumeration::FlavorBase const &_fl )
        { return enumeration::integer_to_structure(_out, _x, _fl); }
                                                        
    void expose_bitset()
    {
      bp::def("get_index", &enumeration::get_index);
      bp::def("count_flavors", &enumeration::count_flavors);
      bp::def("create_flavorbase", &enumeration::create_flavor_base, 
              (bp::arg("card"), bp::arg("nflavors")), "Creates a basis nflavors^k, ( 0<=k<card)." );
      expose_vector<size_t>("FlavorBase", "A basis k^m");
      bp::register_ptr_to_python< boost::shared_ptr< std::vector<size_t> > >();

      bp::def("as_bitstring", &integer_to_bitstring, (bp::arg("integer"), bp::arg("flavorbase")),
              "Converts an integer to a bitstring using flavorbase as  the basis.");
      bp::def("as_numpy", &integer_to_vector, (bp::arg("integer"), bp::arg("flavorbase")),
              "Converts an integer to a numpy array of unsigned integers "
              " using flavorbase as  the basis.");

      bp::def
      (
        "as_structure", 
        &as_structure<Crystal::Structure>,
        ( bp::arg("structure"), bp::arg("x"), bp::arg("flavorbase") ),
        "Fills structure sites using index x. Structure must have correct number of atoms."
      );
      bp::def
      (
        "as_structure", 
        &as_structure< Crystal::TStructure<std::string> >,
        ( bp::arg("structure"), bp::arg("x"), bp::arg("flavorbase") ),
        "Fills structure sites using index x. Structure must have correct number of atoms."
      );
      
      typedef boost::make_signed<enumeration::t_uint>::type t_int;
      bp::scope scope = bp::class_<enumeration::Database>
      (
        "Database", 
        "A large bitset"
      ).def( bp::init<enumeration::Database const&>() )
       .def( "__init__", bp::make_constructor(&create))
       .def( "__init__", bp::make_constructor(&create1<t_int, t_int>))
       .def( "__init__", bp::make_constructor(&create1<enumeration::t_uint, t_int>))
       .def( "__init__", bp::make_constructor(&create1<t_int, enumeration::t_uint>))
       .def("__getitem__", &get_item)
       .def("__setitem__", &set_item)
       .def("__len__", &enumeration::Database::size);

      bp::register_ptr_to_python< boost::shared_ptr<enumeration::Database> >();
    }
  }
} // namespace LaDa

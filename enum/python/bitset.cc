#include "LaDaConfig.h"

#include <sstream>
#include <complex>

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
#include <python/numpy_types.h>
#include <crystal/lattice.h>

#include "bitset.hpp"
#include "../numeric_type.h"
#include "../exceptions.h"



namespace LaDa
{
  
  namespace Python
  {
    namespace bp = boost::python;
    namespace e = LaDa::enumeration;
    bool get_item( e::Database const &_d, e::t_uint _x ) { return _d[_x]; }
    void set_item( e::Database &_d, e::t_uint _x, bool _b ) { _d[_x] = _b; }

    bp::object integer_to_vector(e::t_uint _x, e::FlavorBase const& _fl) 
    {
      std::vector<size_t> result;
      e::integer_to_vector(_x, _fl, result);
      return math::numpy::copy_1darray(result);
    }

    std::string integer_to_bitstring(e::t_uint _x, e::FlavorBase const &_fl )
    {
      try
      {
        return e::integer_to_bitstring(_x, _fl);
      }
      catch(e::integer_too_large &_e)
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

    boost::shared_ptr<e::Database> create(size_t const &_a, size_t const &_b)
    {
      if( _a == 0 or _b == 0 )
      {
        PyErr_SetString(PyExc_ValueError, "Invalid argument: value of 0.\n" );
        bp::throw_error_already_set();
        return boost::shared_ptr<e::Database>();
      }
      try { return e::create_database(_a, _b); }
      catch(e::supercell_too_large &_e)
      {
        std::ostringstream sstr;
        sstr << "Supercell configuration-space too large.\n"
             << "Maximum size is 2^" << sizeof(e::t_uint) << "\n"
             << boost::diagnostic_information(_e);
        PyErr_SetString(PyExc_RuntimeError, sstr.str().c_str());
        bp::throw_error_already_set();
        return boost::shared_ptr<e::Database>();
      }
      catch(boost::exception &_e)
      {
        std::ostringstream sstr;
        sstr << "Error from boost library.\n"
             << boost::diagnostic_information(_e);
        PyErr_SetString(PyExc_RuntimeError, sstr.str().c_str());
        bp::throw_error_already_set();
        return boost::shared_ptr<e::Database>();
      }
      catch(...)
      {
        PyErr_SetString(PyExc_ValueError, "Unknown error while creating database.\n");
        bp::throw_error_already_set();
        return boost::shared_ptr<e::Database>();
      }
    }

    template<class T> typename boost::enable_if<boost::is_signed<T>, bool>::type
      negative( T const &_a ) { return _a < 0; } 
    template<class T> typename boost::disable_if<boost::is_signed<T>, bool>::type
      negative( T const &_a ) { return false; }
    template<class T> typename boost::enable_if<boost::is_signed<T>, e::t_uint>::type
      abs( T const &_a ) { return std::abs(_a); }
    template<class T> typename boost::disable_if<boost::is_signed<T>, e::t_uint>::type
      abs( T const &_a ) { return _a; }

    template< class T1, class T2>
      boost::shared_ptr<e::Database> create1(T1 const &_a, T2 const &_b)
      {
        if( negative(_a) or negative(_b) )
        {
          PyErr_SetString(PyExc_ValueError, "Invalid argument: negative value.\n" );
          bp::throw_error_already_set();
          return boost::shared_ptr<e::Database>();
        }
        return create(_a, _b);
      }


    template<class T_STR> 
      inline void as_structure( Crystal::Lattice & lattice,
                                T_STR &_out, e::t_uint _x,
                                e::FlavorBase const &_fl )
      {
         Crystal::Lattice *old = _out.lattice;
         _out.lattice = &lattice;
         e::integer_to_structure(_out, _x, _fl); 
         _out.lattice = old;
      }

    struct flavor_iter
    {
      flavor_iter   (e::FlavorBase const &_fl, e::t_uint const _x) 
                  : i_it(_fl.rbegin()), i_it_end(_fl.rend()), x(_x), is_first(true) 
      {
        if( x < *i_it * _fl[1] ) return;
        
        PyErr_SetString(PyExc_ValueError, "integer argument is too large.\n");
        boost::python::throw_error_already_set();
        i_it = i_it_end;
        is_first = false;
      }
      flavor_iter   (flavor_iter const &_c) 
                  : i_it(_c.i_it), i_it_end(_c.i_it_end), x(_c.x), is_first(_c.is_first) {}

      flavor_iter iter() { return *this; }
      e::t_uint next()
      {
        if( is_first ) is_first = false;
        else if(i_it != i_it_end) ++i_it;

        if( i_it == i_it_end )
        {
          PyErr_SetString(PyExc_StopIteration, "Stop iteration.\n");
          boost::python::throw_error_already_set();
          return 0u;
        }
        
        e::t_uint const flavor( x / (*i_it) );
        x %= (*i_it);
        return flavor;
      };

      e::t_uint x;
      e::FlavorBase::const_reverse_iterator i_it;
      e::FlavorBase::const_reverse_iterator i_it_end;
      bool is_first;
    };

    boost::shared_ptr<e::FlavorBase> create_flavor_base(size_t _card, size_t _n_flavor)
    {
      try { return e::create_flavor_base(_card, _n_flavor); }
      catch(e::supercell_too_large const &_e)
      {
        std::ostringstream sstr;
        sstr << "Supercell configuration-space too large.\n"
             << "Maximum size is 2^" << sizeof(e::t_uint);
        if( std::string const * mi=boost::get_error_info<e::error_string>(_e) )
          PyErr_SetString( PyExc_ValueError, (sstr.str() + ": " + (*mi)).c_str());
        else
          PyErr_SetString( PyExc_ValueError, (sstr.str() + ".").c_str());
        bp::throw_error_already_set();
        return boost::shared_ptr<e::FlavorBase>();
      }
    }

    void expose_bitset()
    {
      import_array();
      bp::def("get_index", &e::get_index);
      bp::def("count_flavors", &e::count_flavors);
      bp::def("create_flavorbase", &e::create_flavor_base, 
              (bp::arg("card"), bp::arg("nflavors")), "Creates a basis nflavors^k, ( 0<=k<card)." );
      expose_vector<size_t>("FlavorBase", "A basis k^m");
      bp::register_ptr_to_python< boost::shared_ptr< std::vector<size_t> > >();
      
      bp::class_<flavor_iter>("IntegerIterator", "Iterates over atoms of an integer structure.\n",
                              bp::init<e::FlavorBase const &, e::t_uint>())
        .def("__iter__", &flavor_iter::iter )
        .def("next", &flavor_iter::next );

      bp::def("as_bitstring", &integer_to_bitstring, (bp::arg("integer"), bp::arg("flavorbase")),
              "Converts an integer to a bitstring using flavorbase as the basis.");
      bp::def("as_numpy", &integer_to_vector, (bp::arg("integer"), bp::arg("flavorbase")),
              "Converts an integer to a numpy array of unsigned integers "
              " using flavorbase as the basis.");

      bp::def
      (
        "as_structure", 
        &as_structure< Crystal::TStructure<std::string> >,
        (bp::arg("lattice"), bp::arg("structure"), bp::arg("x"), bp::arg("flavorbase") )
      );
      bp::def
      (
        "as_structure", 
        &as_structure<Crystal::Structure>,
        (bp::arg("lattice"), bp::arg("structure"), bp::arg("x"), bp::arg("flavorbase") ),
        "Fills structure sites using index x.\n\n"
        ":Parameters:\n"
        "  lattice : `crystal.Lattice`\n"
        "    Lattice forming the backbone of the enumeration.\n"
        "  structure : `crystal.Structure` or `crystal.rStructure`\n"
        "    a crystal structure with correct number of atoms and lattice vectors.\n"
        "    It will be decorated according to x on output.\n"
        "  x : integer\n    index of the particular decoration.\n"
        "  flavorbase : `_enumeration.FlavorBase`\n" 
        "    the return from `_enumeration.create_flavorbase`.\n"
      );
      
      typedef boost::make_signed<e::t_uint>::type t_int;
      bp::scope scope = bp::class_<e::Database>
      (
        "Database", 
        "A large bitset"
      ).def( bp::init<e::Database const&>() )
       .def( "__init__", bp::make_constructor(&create))
       .def( "__init__", bp::make_constructor(&create1<t_int, t_int>))
       .def( "__init__", bp::make_constructor(&create1<e::t_uint, t_int>))
       .def( "__init__", bp::make_constructor(&create1<t_int, e::t_uint>))
       .def("__getitem__", &get_item)
       .def("__setitem__", &set_item)
       .def("__len__", &e::Database::size);

      bp::register_ptr_to_python< boost::shared_ptr<e::Database> >();
    }
  }
} // namespace LaDa

//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <sstream>

#include <pyublas/numpy.hpp>
#include <boost/python/class.hpp>
#include <boost/python/tuple.hpp>
#include <boost/numeric/ublas/io.hpp>


#include <crystal/structure.h>

#include "../collapse/collapse.h"
#include "../sum_of_separables.h"


namespace LaDa
{
  namespace Python
  {
    namespace bp = boost::python;

    typedef LaDa::atomic_potential::collapse::Collapse Collapse;
    typedef LaDa::atomic_potential::matrix_type matrix_type;
    typedef LaDa::atomic_potential::vector_type vector_type;
    typedef LaDa::atomic_potential::numeric_type numeric_type;

    bp::tuple lsq_data( Collapse const &_coll,
                                   size_t _i) 
    {
      pyublas::numpy_matrix<matrix_type::value_type> mat; 
      pyublas::numpy_vector<vector_type::value_type> vec; 
      _coll.lsq_data( mat, vec, _i );
      return bp::make_tuple( mat, vec );
    }

    numeric_type set_weight(Collapse &_coll, types::t_int _i, numeric_type _w)
    {
      size_t i(_i);
      types::t_int const nbs( _coll.nb_structures() );
      if(_i >= nbs or _i <= -nbs)
      {
        std::ostringstream sstr("Structure index out of range: ");
        sstr << _i << " not int [" << -nbs  << ", " << nbs << "].\n";
        PyErr_SetString(PyExc_IndexError, sstr.str().c_str());
        bp::throw_error_already_set();
        return -1;
      }
      if( _i < 0  ) i = size_t( types::t_int( + _i) );
      return _coll.set_weight(i, _w);
    }

    numeric_type clear_weight(Collapse &_coll, types::t_int _i)
      { return set_weight(_coll, _i, 0e0); }

    pyublas::numpy_vector<vector_type::value_type> coefficients( Collapse const &_coll, size_t _i )
      { return _coll.coefficients(_i); }

    void expose_collapse()
    {
      bp::class_<Collapse>
      ( 
        "Collapse", 
        "Computes matrices for alternating-least-square fit.",
        bp::init<atomic_potential::SumOfSeparables&>() 
      ).def(bp::init<Collapse const&>())
       .def("lsq_data", &lsq_data)
       .def("update", &Collapse::update)
       .def("add", &Collapse::add)
       .def("reassign", &Collapse::reassign)
       .def("coefficients", &coefficients)
       .def("get_weight", &Collapse::get_weight)
       .def("set_weight", &set_weight)
       .def("set_weight", &clear_weight)
       .add_property("values", &Collapse::values_)
       .add_property("scales", &Collapse::scaling_factors_)
       .add_property("fitting_set", &Collapse::fitting_set_)
       .add_property("nb_funcs", &Collapse::nb_funcs_)
       .add_property("nb_coordinates", &Collapse::nb_coordinates)
       .add_property("nb_structures", &Collapse::nb_structures)
       .add_property("y_squared", &Collapse::y_squared)
       .add_property("sum_w", &Collapse::sum_w)
       .add_property("convergence", &Collapse::convergence);
    }

  }
} // namespace LaDa

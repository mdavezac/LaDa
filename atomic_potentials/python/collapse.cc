//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <pyublas/numpy.hpp>
#include <boost/python/class.hpp>
#include <boost/python/tuple.hpp>
#include <crystal/structure.h>

#include "../collapse/collapse.h"
#include "../sum_of_separables.h"


namespace LaDa
{
  namespace Python
  {
    typedef LaDa::atomic_potential::collapse::Collapse Collapse;
    typedef LaDa::atomic_potential::matrix_type matrix_type;
    typedef LaDa::atomic_potential::vector_type vector_type;

    boost::python::tuple lsq_data( Collapse const &_coll,
                                   pyublas::numpy_matrix<matrix_type::value_type> _mat, 
                                   pyublas::numpy_vector<vector_type::value_type> _vec, 
                                   size_t _i) 
    {
      _coll.lsq_data( _mat, _vec, _i );
      return boost::python::make_tuple( _mat, _vec );
    }

    void expose_collapse()
    {
      namespace bp = boost::python;
//     typedef atomic_potential::matrix_type matrix_type;
//     bp::class_<matrix_type>
//     (
//       "Matrix", 
//       "Blind matrix object for the fitting a sum of separable function."
//     ).def( bp::init<matrix_type const&>() );
//     typedef atomic_potential::vector_type vector_type;
//     bp::class_<vector_type>
//     (
//       "Vector", 
//       "Blind vector object for the fitting a sum of separable function."
//     ).def( bp::init<vector_type const&>() );




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
       .add_property("values", &Collapse::values_)
       .add_property("fitting_set", &Collapse::fitting_set_);
    }

  }
} // namespace LaDa

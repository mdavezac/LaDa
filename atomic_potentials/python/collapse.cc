//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

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
    typedef LaDa::atomic_potential::collapse::Collapse Collapse;
    typedef LaDa::atomic_potential::matrix_type matrix_type;
    typedef LaDa::atomic_potential::vector_type vector_type;

//   boost::python::tuple lsq_data( Collapse const &_coll,
//                                  pyublas::numpy_matrix<matrix_type::value_type> _mat, 
//                                  pyublas::numpy_vector<vector_type::value_type> _vec, 
//                                  size_t _i) 
//   {
//     _coll.lsq_data( _mat, _vec, _i );
//     return boost::python::make_tuple( _mat, _vec );
//   }
    boost::python::tuple lsq_data( Collapse const &_coll,
                                   size_t _i) 
    {
      pyublas::numpy_matrix<matrix_type::value_type> mat; 
      pyublas::numpy_vector<vector_type::value_type> vec; 
      _coll.lsq_data( mat, vec, _i );
      return boost::python::make_tuple( mat, vec );
    }

    pyublas::numpy_vector<vector_type::value_type> coefficients( Collapse const &_coll, size_t _i )
      { return _coll.coefficients(_i); }

    void expose_collapse()
    {
      namespace bp = boost::python;

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
       .add_property("values", &Collapse::values_)
       .add_property("scales", &Collapse::scaling_factors_)
       .add_property("fitting_set", &Collapse::fitting_set_)
       .add_property("nb_coordinates", &Collapse::nb_coordinates)
       .add_property("convergence", &Collapse::convergence);
    }

  }
} // namespace LaDa

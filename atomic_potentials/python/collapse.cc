//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python/class.hpp>
#include <crystal/structure.h>

#include "../collapse/collapse.h"
#include "../sum_of_separables.h"


namespace LaDa
{
  namespace Python
  {
    void expose_collapse()
    {
      namespace bp = boost::python;
      typedef atomic_potential::matrix_type matrix_type;
      bp::class_<matrix_type>
      (
        "Matrix", 
        "Blind matrix object for the fitting a sum of separable function."
      ).def( bp::init<matrix_type const&>() );
      typedef atomic_potential::vector_type vector_type;
      bp::class_<vector_type>
      (
        "Vector", 
        "Blind vector object for the fitting a sum of separable function."
      ).def( bp::init<vector_type const&>() );

      typedef LaDa::atomic_potential::collapse::Collapse Collapse;
      bp::class_<Collapse>
      ( 
        "Collapse", 
        "Computes matrices for alternating-least-square fit.",
        bp::init<atomic_potential::SumOfSeparables&>() 
      ).def(bp::init<Collapse const&>())
       .def("lsq_data", &Collapse::lsq_data)
       .def("update", &Collapse::update)
       .def("add", &Collapse::add)
       .def("reassign", &Collapse::reassign)
       .add_property("values", &Collapse::values_)
       .add_property("fitting_set", &Collapse::fitting_set_);
    }

  }
} // namespace LaDa

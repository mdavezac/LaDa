//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python/scope.hpp>
#include <boost/python/class.hpp>

#include "../collapse/fitting_set.h"
#include "../representation.h"
#define LADA_PYTHON_STD_VECTOR_NOPRINT
#include <python/std_vector.hpp>


namespace LaDa
{
  namespace Python
  {
    void expose_fittingset()
    {
      namespace bp = boost::python;
      typedef LaDa::atomic_potential::collapse::FittingSet FittingSet;
      
      bp::scope scope = bp::class_<FittingSet>
              ("FittingSet", "Holds fitting related info. For Debugging only.")
        .add_property("coordinates", &FittingSet::coordinates_)
        .add_property("energies", &FittingSet::energies_)
        .add_property("weights", &FittingSet::weights_)
        .def("add", &FittingSet::add);

      expose_vector<FittingSet::t_Coordinates::value_type>
          ("coord_vec0", "Implementation only.");
      expose_vector<FittingSet::t_Coordinates::value_type::value_type>
          ("coord_vec1", "Implementation only.");
      expose_vector<FittingSet::t_Coordinates::value_type::value_type::value_type>
          ("coord_vec2", "Implementation only.");
      
      bp::class_< FittingSet::t_Weights::value_type >("weight_pair")
        .add_property("fitting", &FittingSet::t_Weights::value_type::first)
        .add_property("representation", &FittingSet::t_Weights::value_type::second);
 
      expose_vector<FittingSet::t_Weights::value_type>("vec_weights", "Implementation only.");
 
    }

  }
} // namespace LaDa

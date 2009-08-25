//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python/class.hpp>

#include "../collapse/fitting_set.h"
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
      
      expose_vector<FittingSet::t_AtomType>("FittingSet_AtomType", "Implementation only.");
      expose_vector<std::vector<FittingSet::t_Energies> >("FittingSet_RepCoord",
                                                          "Implementation only.");
      expose_vector<std::vector<std::vector<FittingSet::t_Energies> > >
            ("FittingSet_StrCoord", "Implementation only.");
      expose_vector<std::vector<std::vector<std::vector<FittingSet::t_Energies> > > >
            ("FittingSet_Coord", "Implementation only.");
 
      expose_vector<FittingSet::t_Energies>("FittingSet_NumericVector", "Implementation only.");
 
      bp::class_< FittingSet::t_Weights::value_type >("FittingSet_WeightsPair")
        .add_property("str_weight", &FittingSet::t_Weights::value_type::first)
        .add_property("rep_weights", &FittingSet::t_Weights::value_type::second);
 
      expose_vector<FittingSet::t_Weights>("FittingSet_Weights", "Implementation only.");
 
      bp::class_<FittingSet>("FittingSet", "Holds fitting related info. For Debugging only.")
        .add_property("coordinates", &FittingSet::coordinates_)
        .add_property("energies", &FittingSet::energies_)
        .add_property("weights", &FittingSet::weights_);
    }

  }
} // namespace LaDa

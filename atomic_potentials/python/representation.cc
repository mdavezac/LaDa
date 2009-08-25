//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <sstream>
#include <complex>

#include <boost/python/class.hpp>

#include <crystal/structure.h>

#define LADA_PYTHON_STD_VECTOR_NOPRINT
#include <python/std_vector.hpp>
#include <python/misc.hpp>

#include "../representation.h"


namespace LaDa
{
  namespace Python
  {
    void expose_representation()
    {
      namespace bp = boost::python;
      typedef LaDa::atomic_potential::VariableSet::t_Variable Variable;
      bp::class_<Variable>
      ( 
        "Variable", 
        "Coordinate and specie variable from a representations.",
        bp::init<Variable const&>() 
      ).add_property("coordinate", &Variable::first)
       .add_property("specie", &Variable::second);
      
      typedef LaDa::atomic_potential::VariableSet::t_Variables Variables;
      expose_vector<Variables>( "Variables", "Container of variables for a representations.");
        
      typedef LaDa::atomic_potential::VariableSet VariableSet;
      bp::class_<VariableSet>
      ( 
        "VariableSet",
        "Set of variables + weight. A set of VariableSets is a representation.",
        bp::init<VariableSet const&>()
      ).add_property("weight", &VariableSet::weight)
       .add_property("variables", &VariableSet::variables);

      typedef LaDa::atomic_potential::Representation Representation;
      bp::class_<Representation>
      ( 
        "Representation", 
        "Symmetrized sets of variables forming a representation within "
        "the sum of separable functions framework.",
        bp::init<Representation const&>() 
      ).def( bp::init< Crystal::TStructure<std::string> const&, size_t>() )
       .def( "__str__", &tostream<Representation const&> )
       .def( "__len__", &Representation::size)
       .def( "nb_coords", &Representation::nb_coords)
       .def( "nb_atoms", &Representation::nb_atoms);
    }

  }
} // namespace LaDa

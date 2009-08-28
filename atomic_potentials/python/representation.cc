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
    struct RepresentationIter 
    {
      RepresentationIter   ( atomic_potential::Representation const &_rep )
                         : first_(true), cit_(_rep.begin()), cit_end_(_rep.end()) {}
      RepresentationIter   ( RepresentationIter const &_c )
                         : cit_(_c.cit_), cit_end_(_c.cit_end_), first_(_c.first_) {}

      RepresentationIter &iter()  { return *this; }
      atomic_potential::Representation::const_iterator::reference next()
      {
        namespace bp = boost::python;
        if( first_ ) first_ = false; 
        else 
        {
          ++cit_;
          if( cit_ == cit_end_ )
          {
            PyErr_SetString
            (
              PyExc_StopIteration, 
              "Error while computing transform to smith normal form.\n" 
            );
            bp::throw_error_already_set();
            --cit_;
          }
        }
        return *cit_;
      }

      atomic_potential::Representation::const_iterator cit_;
      atomic_potential::Representation::const_iterator cit_end_;
      bool first_;
    };

    RepresentationIter create_iter( atomic_potential::Representation const &_rep )
      { return RepresentationIter(_rep); }

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
      
      typedef LaDa::atomic_potential::VariableSet::t_Variables t_Variables;
      expose_vector<t_Variables::value_type>
        ( "Variables", "Container of variables for a representations.");
        
      typedef LaDa::atomic_potential::VariableSet VariableSet;
      bp::class_<VariableSet>
      ( 
        "VariableSet",
        "Set of variables + weight. A set of VariableSets is a representation.",
        bp::init<VariableSet const&>()
      ).add_property("weight", &VariableSet::weight)
       .add_property("variables", &VariableSet::variables)
       .def("__str__", &tostream<VariableSet>);

      typedef LaDa::atomic_potential::Representation Representation;
      bp::class_<RepresentationIter>
      (
        "RepresentationIterator", 
        "Representation iterator.",
        bp::init<RepresentationIter const&>()
      ).def("__iter__", &RepresentationIter::iter, bp::return_internal_reference<1>() )
       .def("next", &RepresentationIter::next, bp::return_internal_reference<1>() );

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
       .def( "nb_atoms", &Representation::nb_atoms)
       .def( "__iter__", &create_iter);
    }

  }
} // namespace LaDa

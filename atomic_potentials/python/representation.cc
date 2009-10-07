//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <sstream>
#include <complex>

#include <boost/python/class.hpp>
#include <boost/python/str.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/scope.hpp>

#include <crystal/structure.h>

#include <python/std_vector.hpp>
#include <python/misc.hpp>

#include "../representation.h"


namespace LaDa
{
  namespace Python
  {
    typedef LaDa::atomic_potential::Representation Representation;
    size_t getcs () { return Representation::coord_system; }
    struct RepresentationIter 
    {
      RepresentationIter   ( Representation const &_rep )
                         : first_(true), cit_(_rep.begin()), cit_end_(_rep.end()) {}
      RepresentationIter   ( RepresentationIter const &_c )
                         : cit_(_c.cit_), cit_end_(_c.cit_end_), first_(_c.first_) {}

      RepresentationIter &iter()  { return *this; }
      Representation::const_iterator::reference next()
      {
        namespace bp = boost::python;
        if( first_ ) first_ = false; 
        else 
        {
          ++cit_;
          if( cit_ == cit_end_ )
          {
            PyErr_SetString(PyExc_StopIteration, "Iterator out-of-range");
            bp::throw_error_already_set();
            --cit_;
          }
        }
        return *cit_;
      }

      Representation::const_iterator cit_;
      Representation::const_iterator cit_end_;
      bool first_;
    };

    RepresentationIter create_iter( Representation const &_rep )
      { return RepresentationIter(_rep); }

    Representation::const_iterator::reference getitem( Representation const& _rep, size_t _i)
    {
      if( _i >= _rep.size() )
      {
        PyErr_SetString(PyExc_IndexError, "Index out-of-range.");
        boost::python::throw_error_already_set();
        {
          static Representation::const_iterator::value_type val;
          return val;
        }
      }
      return _rep[_i];
    }

    std::ostream& operator<<( std::ostream& _stream,
                              atomic_potential::VariableSet::t_Variable const &_var )
      { return atomic_potential::operator<<(_stream, _var); }

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
       .add_property("specie", &Variable::second)
       .def("__str__", &tostream<atomic_potential::VariableSet::t_Variable>);
      
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

      bp::class_<RepresentationIter>
      (
        "RepresentationIterator", 
        "Representation iterator.",
        bp::init<RepresentationIter const&>()
      ).def("__iter__", &RepresentationIter::iter, bp::return_internal_reference<1>() )
       .def("next", &RepresentationIter::next, bp::return_internal_reference<1>() );

      bp::object rep = bp::class_<Representation>
      ( 
        "Representation", 
        "Symmetrized sets of variables forming a representation within "
        "the sum of separable functions framework.",
        bp::init<Representation const&>() 
      ).def( bp::init< Crystal::TStructure<std::string> const&, size_t>() )
       .def( "__str__", &tostream<Representation> )
       .def( "__len__", &Representation::size)
       .add_property( "nb_coordinates", &Representation::nb_coords)
       .add_property( "nb_atoms", &Representation::nb_atoms)
       .def( "__iter__", &create_iter)
       .def( "__getitem__", &getitem, bp::return_internal_reference<1>());

      bp::scope scope = rep;

      bp::enum_< Representation::CoordinateSystem >( "CoordinateSystem" )
        .value( "cartesian", Representation::cartesian )
        .value( "spherical", Representation::spherical )
        .export_values();
      
      rep.attr("coord_system") =  Representation::coord_system;
    }

  }
} // namespace LaDa

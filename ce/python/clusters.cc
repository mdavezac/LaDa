//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <sstream>
#include <complex>

#include <boost/python/scope.hpp>
#include <boost/python/class.hpp>
#include <boost/python/errors.hpp>
#include <boost/python/register_ptr_to_python.hpp>


#include "cluster.h"


namespace LaDa
{
  namespace Python
  {
    CE::Cluster apply_symmetry( CE::Cluster const &_self, Crystal::SymmetryOperator const &_op );
    {
      CE::Cluster result(_self);
      result.apply_symmetry(_op);
      return result;
    }

    boost::shared_ptr< std::vector<CE::Clusters> > 
      equiv_clusters( CE::Cluster &_self, Crystal::Lattice &_lat )
      {
        boost::shared_ptr< std::vector<CE::Clusters> >
          result( new std::vector<CE::Clusters>(1, _self) );
        add_equivalent_clusters( _lat, *result);
        return result;
      }

    void expose_clusters()
    {
      namespace bp = boost::python;

      bp::scope scope = bp::class_<CE::Cluster>("Cluster", "A cluster (figure).")
          .def(bp::init<Cluster const&>())
          .add_property("eci", &CE::Cluster::eci, "Interaction energy.")
          .add_property("vectors", &Neighbor::vectors, "Spin positions.")
          .def("apply_symmetry", &apply_symmetry, "Returns a transformed cluster.\n")
          .def("__len__", &CE::Cluster::size, "Returns the order of the cluster.\n")
          .def("equivalents", &equivalents, "Returns array of arrays of equivalent clusters.\n");

      expose_vector<CE::Cluster>
        ("Clusters", "An array of (presumably) equivalent clusters.\n");
      bp::register_ptr_to_python
      < 
        boost::shared_ptr< std::vector<CE::Cluster> >
      >();
      expose_vector< std::vector<CE::Cluster> >
        ("Clusters", "An array of arrays of (presumably) equivalent clusters.\n");
    }

  }
} // namespace LaDa

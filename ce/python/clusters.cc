//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <iostream>
#include <sstream>
#include <complex>

#include <boost/python/def.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/class.hpp>
#include <boost/python/errors.hpp>
#include <boost/python/register_ptr_to_python.hpp>

#include "clusters.hpp"
#include "../cluster.h"

#include <python/misc.hpp>
#include <python/std_vector.hpp>


namespace LaDa
{
  namespace Python
  {
    CE::Cluster apply_symmetry( CE::Cluster const &_self, Crystal::SymmetryOperator const &_op )
    {
      CE::Cluster result(_self);
      result.apply_symmetry(_op);
      return result;
    }

    boost::shared_ptr< std::vector<CE::Cluster> > 
      equiv_clusters( CE::Cluster &_self, Crystal::Lattice &_lat )
      {
        boost::shared_ptr< std::vector<CE::Cluster> >
          result( new std::vector<CE::Cluster>(1, _self) );
        add_equivalent_clusters( _lat, *result);
        return result;
      }

    std::ostream &operator<<( std::ostream &_stream, std::vector<CE::Cluster> const &_cls)
    {
      _stream << "Cluster Class:\n";
      foreach( CE::Cluster const &cluster, _cls)
        _stream << cluster;
      return _stream;
    }
    

    struct Rot
    {
      Rot( Crystal::SymmetryOperator const &_op ) : op(_op) {}
      Rot( Rot const &_op ) : op(_op.op) {} 
      Crystal::SymmetryOperator const & op;
      void operator()( atat::rVector3d &_vec ) const { _vec = op.op * _vec; }
      void operator()( CE::Cluster &_cls ) const
        { std::for_each(_cls.vectors.begin(), _cls.vectors.end(), *this); }
      void operator()( CE::t_Clusters &_cls ) const
        { std::for_each(_cls.begin(), _cls.end(), *this); }
      void operator()( CE::t_ClusterClasses &_cls ) const
        { std::for_each(_cls.begin(), _cls.end(), *this); }
    };
    struct Op
    {
      Op( Crystal::SymmetryOperator const &_op ) : op(_op) {}
      Op( Op const &_op ) : op(_op.op) {} 
      Crystal::SymmetryOperator const & op;
      void operator()( CE::Cluster &_cls ) const { _cls.apply_symmetry(op); }
      void operator()( CE::t_Clusters &_cls ) const
        { std::for_each(_cls.begin(), _cls.end(), *this); }
      void operator()( CE::t_ClusterClasses &_cls ) const
        { std::for_each(_cls.begin(), _cls.end(), *this); }
    };

    template< class T >
      boost::shared_ptr<T> apply_sym(T const &_self, Crystal::SymmetryOperator const &_op)
      {
        boost::shared_ptr<T> result( new T(_self) );
        std::for_each(result->begin(), result->end(), Op(_op));
        return result;
      }
    template< class T >
      boost::shared_ptr<T> apply_rot(T const &_self, Crystal::SymmetryOperator const &_op)
      {
        boost::shared_ptr<T> result( new T(_self) );
        std::for_each(result->begin(), result->end(), Rot(_op));
        return result;
      }
    CE::Cluster apply_rot0(CE::Cluster const &_self, Crystal::SymmetryOperator const &_op)
    {
      CE::Cluster result(_self);
      Rot o(_op); o(result);
      return result;
    }


    void expose_clusters()
    {
      namespace bp = boost::python;
      bp::def("apply_rotation", &apply_rot0, "Apply rotation operation to a cluster.\n");
      bp::def("apply_rotation", &apply_rot<CE::t_Clusters>,
              "Apply rotation operation to an array of equivalent clusters.\n");
      bp::def("apply_rotation", &apply_rot<CE::t_ClusterClasses>,
              "Apply rotation operation to an array of classes of equivalent clusters.\n");

      bp::def("apply_symmetry", &apply_symmetry, "Apply symmetry operation to single cluster.\n");
      bp::def("apply_symmetry", &apply_sym<CE::t_Clusters>,
              "Apply symmetry operation to an array of equivalent clusters.\n");
      bp::def("apply_symmetry", &apply_sym<CE::t_ClusterClasses>,
              "Apply symmetry operation to an array of classes of equivalent clusters.\n");

      bp::scope scope = bp::class_<CE::Cluster>("Cluster", "A cluster (figure).")
          .def(bp::init<CE::Cluster const&>())
          .add_property("eci", &CE::Cluster::eci, "Interaction energy.")
          .add_property("vectors", &CE::Cluster::vectors, "Spin positions.")
          .def("apply_symmetry", &apply_symmetry, "Returns a transformed cluster.\n")
          .def("apply_rotation", &apply_rot0, "Returns a rotated cluster.\n")
          .def("__len__", &CE::Cluster::size, "Returns the order of the cluster.\n")
          .def("__str__", &tostream<CE::Cluster>, "Dumps the cluster to a string.\n")
          .def("equivalents", &equiv_clusters, "Returns array of arrays of equivalent clusters.\n");


      expose_vector<CE::Cluster>
        ("Clusters", "An array of (presumably) equivalent clusters.\n")
        .def("apply_rotation", &apply_rot<CE::t_Clusters>,
             "Returns array of rotated clusters.")
        .def("apply_symmetry", &apply_sym<CE::t_Clusters>,
             "Returns array of transformed clusters.");

      bp::register_ptr_to_python
      < 
        boost::shared_ptr< std::vector<CE::Cluster> >
      >();
 
      expose_vector< std::vector<CE::Cluster> >
        ("Clusters", "An array of arrays of (presumably) equivalent clusters.\n")
        .def("apply_rotation", &apply_rot<CE::t_ClusterClasses>,
             "Returns array of rotated clusters.")
        .def("apply_symmetry", &apply_sym<CE::t_ClusterClasses>,
             "Returns array of transformed clusters.");

      bp::register_ptr_to_python
      < 
        boost::shared_ptr< std::vector< std::vector<CE::Cluster> > >
      >();
    }

  }
} // namespace LaDa

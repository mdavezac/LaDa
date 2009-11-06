//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<boost/lexical_cast.hpp>

#include <print/manip.h>

#include "mlclusters.h"


namespace LaDa
{
  namespace CE 
  {
    void MLClusters :: init(MLCluster const &_cls, types::t_real _eci)
    {
      eci = _eci;
      resize(1);
      back() = _cls;

      // checks for existence of lattice.
      LADA_DOASSERT(Crystal::Structure::lattice,"No lattice found.\n");
      Crystal::Lattice lattice = *Crystal::Structure::lattice;
      if( lattice.space_group.size() == 0 ) lattice.find_space_group();
      if( lattice.space_group.size() == 0 ) return;

      Crystal::Lattice::t_SpaceGroup :: const_iterator i_op = lattice.space_group.begin();
      Crystal::Lattice::t_SpaceGroup :: const_iterator const i_op_end = lattice.space_group.end();
      for(; i_op != i_op_end; ++i_op)
      {
        MLCluster cluster = _cls; 
        cluster.apply_symmetry(*i_op);
        if( end() == std::find(begin(), end(), cluster) ) push_back(cluster);
      }
    }

    bool MLClusters :: Load( const TiXmlElement &_node )
    {
      clear();
      if( not _node.Attribute("eci") ) eci = 0e0;
      else eci = boost::lexical_cast<types::t_real>(_node.Attribute("eci"));

      const TiXmlElement *child = _node.FirstChildElement( "Cluster" );
      MLCluster cluster;
      for ( ; child; child=child->NextSiblingElement( "Cluster" ) )
      {
        if( not cluster.Load(*child) ) return false;
        push_back(cluster);
      }

      return true;
    }

    bool MLClusters::operator==(MLClusters const &_b) const
    {
      if( size() != _b.size() ) return false;
      MLClusters::const_iterator i_cls = _b.begin();
      MLClusters::const_iterator const i_cls_end = _b.end();
      for(; i_cls != i_cls_end; ++i_cls)
        if( end() == std::find(begin(), end(), *i_cls) ) return false;
      return true;
    }
    
    types::t_real MLClusters::operator()( Crystal::Structure const &_str, 
                                          std::vector< std::vector<size_t> > const &_map,
                                          Crystal::t_SmithTransform const &_transform ) const
    {
      if( empty() ) return eci;

      types::t_real result(0);
      const_iterator i_cluster = begin();
      const_iterator const i_cluster_end = end();
      for(; i_cluster != i_cluster_end; ++i_cluster) // loop over clusters classes.
        result += (*i_cluster)(_str, _map, _transform);
      return result;
    }

    std::ostream& operator<<( std::ostream &_stream,  MLClusters const &_class )
    {
      if( _class.size() == 0 )
        return _stream << "Cluster class J0: " << _class.eci << "\n";

      _stream << "Cluster Class: " << _class.eci << "\n";
      std::ostringstream stream;
      MLClusters::const_iterator i_cls = _class.begin();
      MLClusters::const_iterator const i_cls_end = _class.end();
      for(; i_cls != i_cls_end; ++i_cls) stream << *i_cls;
      return _stream << Print::indent( stream.str(), "  ");
    }
    
  } // namespace CE
}

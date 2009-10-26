//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <algorithm>

#include <crystal/which_site.h>
#include <crystal/neighbors.h>
#include <opt/ndim_iterator.h>

#include "create_pairs.h"
#include "create_clusters.h"

namespace LaDa
{
  namespace CE
  {

    //! \brief Creates all clusters of oreder \a _order up to \a _n neighbors.
    boost::shared_ptr<t_ClusterClasses>
      create_clusters( Crystal::Lattice const &_lat, size_t _order,
                       size_t _neighbor, size_t _origin )
      {
        boost::shared_ptr<t_ClusterClasses> result(new t_ClusterClasses);
        create_clusters(*result, _lat, _order, _neighbor, _origin);
        return result;
      }

    //! \brief Creates all clusters of oreder \a _order up to \a _n neighbors.
    void create_clusters( t_ClusterClasses &_out, Crystal::Lattice const &_lat,
                          size_t _order, size_t _neighbor, size_t _origin )
    {
      if( _order < 2 ) 
      {
        CE::Cluster cls; 
        cls.eci = 0e0;
        _out.push_back( t_Clusters(1, cls) );
        if( _order == 1 ) cls.vectors.push_back( atat::rVector3d(0,0,0) );
        return;
      };

      LADA_ASSERT(_origin < _lat.sites.size(), "Site index out of range.\n");
      atat::rVector3d const &origin = _lat.sites[_origin].pos;


      if( _neighbor == 0 ) return;

      // Computes point group at lattice point.
      std::vector<Crystal::SymmetryOperator> point_group;
      foreach(Crystal::SymmetryOperator const &op, _lat.space_group)
        if( Crystal::which_site(op(origin), !_lat.cell, _lat.sites) == _origin )
          point_group.push_back(op);
          
      // Creates a typical cluster.
      Cluster cluster;
      cluster.vectors.resize(_order);
      cluster.vectors[0] = atat::rVector3d(0,0,0);
      cluster.eci = 0e0;

      // gets all neighbors.
      size_t const N(_out.size()+_neighbor);
      size_t size_try;
      switch( _neighbor ) // includes enough neighbors for fcc or bcc. Hence for any lattice?
      {
        default: size_try = 30*_neighbor; break;
        case  1: size_try =  12; break;
        case  2: size_try =  18; break;
        case  3: size_try =  42; break;
        case  4: size_try =  54; break;
        case  5: size_try =  78; break;
        case  6: size_try =  88; break;
        case  7: size_try = 134; break;
        case  8: size_try = 140; break;
        case  9: size_try = 164; break;
        case 10: size_try = 176; break;
        case 11: size_try = 200; break;
      };
      
      std::vector<atat::rVector3d> neighbors;
      do
      { // creates vector of neighbors
        Crystal::Neighbors neighbors_(size_try, origin);
        Crystal::Neighbors::const_iterator i_first = neighbors_.begin(_lat);
        Crystal::Neighbors::const_iterator i_begin = i_first;
        Crystal::Neighbors::const_iterator const i_end = neighbors_.end();
        size_t s(0), n(0);
        types::t_real current_norm(0);
        for(; i_first != i_end; ++i_first, ++n)
        {
          if( Fuzzy::is_zero(current_norm - i_first->distance) ) continue;
          if( s == _neighbor ) break;
          current_norm = i_first->distance;
          ++s;
        }
        if( s != _neighbor ) continue;
        neighbors.reserve(n);
        for(; i_begin != i_first; ++i_begin )
          neighbors.push_back(i_begin->pos);
        break;
      }
      while(true);
      // creates an order-dimensional iterator 
      opt::NDimIterator< size_t, std::less<size_t> > iterator;
      for( size_t i(1); i < _order; ++i )
        iterator.add( 0, neighbors.size() );

      do
      {
        // performs loop over i < j < k < ... only to avoid double counting.
        bool do_iterate = false;
        for(size_t i(1); i < _order; ++i)
        {
          if( iterator.access(i) > iterator.access(i-1) ) continue;
          
          do_iterate = true;
          break;
        }
        if( do_iterate ) continue;

        // creates cluster.
        for(size_t i(1); i < _order; ++i)
          cluster.vectors[i] = neighbors[iterator.access(i-1)]; 

        // now checks whether it should be part of a cluster class somewhere.
        std::vector<Crystal::SymmetryOperator> :: const_iterator i_op = point_group.begin();       
        std::vector<Crystal::SymmetryOperator> :: const_iterator const i_op_end = point_group.end();       
        t_ClusterClasses::iterator i_class = _out.begin();
        t_ClusterClasses::iterator const i_class_end = _out.end();
        for(; i_op != i_op_end; ++i_op );
        {
          Cluster transformed( cluster ); transformed.apply_symmetry( *i_op ); 

          for(i_class = _out.begin(); i_class != i_class_end; ++i_class )
          {
            if( i_class->size() == 0 ) break;
            if( i_class->front().size() != _order ) continue;

            if( i_class->end() != std::find(i_class->begin(), i_class->end(), cluster) ) 
              break;
          }

          // if cluster found, break.
          if( i_class != i_class_end ) break;
        }

        // adds new cluster.
        if( i_op == i_op_end )  _out.push_back( t_Clusters(1, cluster) );
        else i_class->push_back(cluster);

      } while(++iterator);
    }
  } // end of namespace CE
} // namespace LaDa

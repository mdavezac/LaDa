//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <algorithm>
#include <iterator>

#include <crystal/which_site.h>
#include <crystal/neighbors.h>
#include <opt/ndim_iterator.h>

#include "create_pairs.h"
#include "create_clusters.h"

namespace LaDa
{
  namespace CE
  {
    boost::shared_ptr<t_MLClusterClasses>
      create_clusters( Crystal::Lattice const &_lat, size_t const _order,
                       size_t const _neighbor, size_t const _origin )
      {
        boost::shared_ptr<t_MLClusterClasses> result(new t_MLClusterClasses);
        create_clusters(*result, _lat, _order, _neighbor, _origin);
        return result;
      }

    void create_clusters( t_MLClusterClasses &_out, Crystal::Lattice const &_lat,
                          size_t const _order, size_t const _neighbor, size_t const _origin )
    {
      _out.clear();

      // J0
      MLClusters clusters; 
      clusters.eci = 0e0;
      if( _order == 0u )
      {
        _out.push_back( clusters );
        return;
      }
      
      LADA_ASSERT(_origin < _lat.sites.size(), "Site index out of range.\n");
      // Creates a typical cluster.
      CE::MLCluster cluster; 
      cluster.origin.site = _origin;
      cluster.origin.pos = _lat.sites[_origin].pos;
      if( _order == 1u )  // J1
      {
        clusters.init(cluster);
        _out.push_back(clusters);
        return;
      };

      LADA_ASSERT( _neighbor != 0u, "Incoherent arguments.\n");
      cluster.resize(_order-1);

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
      
      std::vector<Crystal::Neighbor> neighbors;
      do
      { // creates vector of neighbors
        std::cout << "size_try: " << size_try << "\n";
        Crystal::Neighbors neighbors_(size_try, cluster.origin.pos);
        Crystal::Neighbors::const_iterator i_first = neighbors_.begin(_lat);
        Crystal::Neighbors::const_iterator i_begin = i_first;
        Crystal::Neighbors::const_iterator const i_end = neighbors_.end();
        size_t s(0), n(0);
        types::t_real current_norm(0);
        for(; i_first != i_end; ++i_first, ++n)
        {
          if( _lat.sites[i_first->index].type.size() < 2 ) continue;
          if( Fuzzy::is_zero(current_norm - i_first->distance) ) continue;
          if( s == _neighbor ) break;
          current_norm = i_first->distance;
          ++s;
        }
        if( s != _neighbor ) { size_try <<= 1; continue; } // did not find all neighbors.

        neighbors.reserve(n);
        std::copy(i_begin, i_first, std::back_inserter(neighbors));
        break;
      }
      while(true);
      std::cout << "done neighbors.\n";
      // creates an order-dimensional iterator 
      opt::NDimIterator< size_t, std::less<size_t> > iterator;
      for( size_t i(1); i < _order; ++i )
        iterator.add( 0, neighbors.size() );

      do
      {
        // performs loop over i < j < k < ... only to avoid double counting.
        bool dont_iterate = true;
        for(size_t i(1); i < _order-1 and dont_iterate; ++i)
          if( iterator.access(i) <= iterator.access(i-1) ) dont_iterate = false;
        if( not dont_iterate ) continue;
        std::cout << 1u << "\n";

        // creates cluster.
        for(size_t i(0); i < _order-1; ++i)
        {
          cluster[i].pos = neighbors[iterator.access(i)].pos; 
          cluster[i].site = neighbors[iterator.access(i)].index; 
        }

        std::cout << 2u << "\n";
        // now checks whether it should be part of a cluster class somewhere.
        t_MLClusterClasses::iterator i_class = _out.begin();
        t_MLClusterClasses::iterator const i_class_end = _out.end();
        for(; i_class != i_class_end; ++i_class );
          if( i_class->end() != std::find(i_class->begin(), i_class->end(), cluster) )
            break;
        std::cout << 3u << "\n";

        // adds new cluster.
        if( i_class == i_class_end ) 
        {
        std::cout << 4u << "\n";
          clusters.init( cluster );
          _out.push_back( clusters );
          std::cout << "_out size: " << _out.size() << "\n";
        }
        std::cout << 5u << "\n";

      } while(++iterator);
    }

  } // end of namespace CE
} // namespace LaDa

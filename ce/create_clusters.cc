//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <algorithm>
#include <iterator>
#include <set>

#include <boost/dynamic_bitset.hpp>

#include <crystal/which_site.h>
#include <crystal/neighbors.h>
#include <opt/ndim_iterator.h>
#include <opt/fuzzy.h>
#include <opt/pow.h>

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

    typedef boost::dynamic_bitset<long long unsigned int> Database;
    typedef Database :: size_type t_uint;
    typedef std::vector<Crystal::Neighbor> t_Neighbors;
    typedef std::vector<t_uint> t_Flavors;

    void integer_to_cluster( t_uint _x, 
                             t_Flavors const &_flavors, 
                             t_Neighbors const &_neighbors,
                             MLCluster &_out )
    {
      _out.clear();
      t_Flavors :: const_reverse_iterator i_flavor = _flavors.rbegin();
      t_Flavors :: const_reverse_iterator i_flavor_end = _flavors.rend();
      for(; i_flavor != i_flavor_end; ++i_flavor)
      {
         LADA_ASSERT(*i_flavor != 0, "Internal bug.\n");
         t_uint const flavor( _x / (*i_flavor) );
         _x %= (*i_flavor);
         LADA_ASSERT(flavor < _neighbors.size(), "Index out of range.\n");
         t_Neighbors::value_type const &neighbor( _neighbors[flavor] );
         MLCluster::Spin const spin = {neighbor.index, neighbor.pos};
         _out.push_back(spin);
      }
    }

    t_uint cluster_to_integer( MLCluster const &_in,
                               t_Flavors const &_flavors, 
                               t_Neighbors const &_neighbors )
    {
      LADA_ASSERT(_flavors.size() == _in.size(), "Flavors and cluster have different sizes.\n")
      MLCluster::const_iterator i_vec = _in.begin(); 
      MLCluster::const_iterator const i_vec_end = _in.end(); 
      std::set<t_uint> numbers;
      for(; i_vec != i_vec_end; ++i_vec)
      {
        t_Neighbors :: const_iterator i_neigh = _neighbors.begin();
        t_Neighbors :: const_iterator const i_neigh_end = _neighbors.end();
        t_uint i(0);
        for(; i_neigh != i_neigh_end; ++i_neigh, ++i)
        {
          if( not Fuzzy::is_zero(i_neigh->pos(0) - i_vec->pos(0)) ) continue;
          if( not Fuzzy::is_zero(i_neigh->pos(1) - i_vec->pos(1)) ) continue;
          if( Fuzzy::is_zero(i_neigh->pos(2) - i_vec->pos(2)) ) break;
        }
        // returns "false" if could not find vector.
        if( i_neigh == i_neigh_end )
          return _flavors.back()*_neighbors.size()*_neighbors.size();

        numbers.insert(i);
      }
      t_uint result(0);
      std::set<t_uint>::const_iterator i_n = numbers.begin();
      std::set<t_uint>::const_iterator i_n_end = numbers.end();
      t_Flavors :: const_reverse_iterator i_flavor = _flavors.rbegin();
      for(; i_n != i_n_end; ++i_n, ++i_flavor)
        result += (*i_flavor) * (*i_n);
      return result;
    }
    bool permissible( t_uint _x, t_Flavors const &_flavors ) 
    {
      t_Flavors :: const_reverse_iterator i_flavor = _flavors.rbegin();
      t_Flavors :: const_reverse_iterator const i_flavor_end = _flavors.rend();
      std::set<t_uint> numbers;
      for(; i_flavor != i_flavor_end; ++i_flavor)
      {
        t_uint const flavor( _x / (*i_flavor) );
        _x %= (*i_flavor);
        if( numbers.count(flavor) ) return false;
        numbers.insert(flavor);
        if( *numbers.rbegin() != flavor ) return false;
      }
      return true;
    }

    void create_flavors(t_uint _order, t_uint const _base, t_Flavors &_out)
    {
      _out.clear();
      for(t_uint i(0), j(1); i < _order - 1; ++i, j*=_base)
        _out.push_back(j);
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
        case  1: size_try =  13; break;
        case  2: size_try =  19; break;
        case  3: size_try =  43; break;
        case  4: size_try =  55; break;
        case  5: size_try =  79; break;
        case  6: size_try =  89; break;
        case  7: size_try = 135; break;
        case  8: size_try = 141; break;
        case  9: size_try = 165; break;
        case 10: size_try = 177; break;
        case 11: size_try = 201; break;
      };
      
      t_Neighbors neighbors;
      do
      { // creates vector of neighbors
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
        if( i_first == i_end ) { size_try <<= 1; continue; } // did not find all neighbors.

        neighbors.reserve(n);
        std::copy(i_begin, i_first, std::back_inserter(neighbors));
        break;
      }
      while(true);

      // creates flavors and database.
      t_uint const Nmax(opt::pow(neighbors.size(), _order-1));
      Database database(Nmax);
      database.flip();
      t_Flavors flavors;
      create_flavors(_order, neighbors.size(), flavors);

      // database refers to all possible clusters.
      for(t_uint n(0); n < Nmax; ++n)
      {
        // checks that this cluster is not already known.
        if( not database[n] ) continue;
        // checks if permissible number.
        if( not permissible(n, flavors) ) continue;

        // initializes cluster.
        integer_to_cluster(n, flavors, neighbors, cluster);

        // Initializes cluster
        clusters.init( cluster );

        // now turns off these clusters in the database,
        // and makes sure cluster fits in shell.
        // Can we do this when determining neighbors? 
        // Not sure for multi-site lattices.
        bool is_in_shell(true);
        MLClusters :: const_iterator i_cluster = clusters.begin();
        MLClusters :: const_iterator const i_cluster_end = clusters.end();
        for(; i_cluster != i_cluster_end; ++i_cluster)
        {
          // original site.
          if(i_cluster->origin.site == _origin) 
          {
            t_uint const n_prime = cluster_to_integer(*i_cluster, flavors, neighbors);
            if( n_prime < Nmax ) database[n_prime] = false;
            else is_in_shell = false;
          }
          // shifted sites.
          MLCluster::t_Spins :: const_iterator i_spin = i_cluster->begin();
          MLCluster::t_Spins :: const_iterator const i_spin_end = i_cluster->end();
          for(; i_spin != i_spin_end; ++i_spin)
          {
            if( _origin != i_spin->site ) continue;
            MLCluster shifted;
            shifted.origin.site = i_spin->site; 
            shifted.origin.pos = _lat.sites[i_spin->site].pos;
            MLCluster::t_Spins :: const_iterator i_spin_add = i_cluster->begin();
            for(; i_spin_add != i_spin_end; ++i_spin_add)
              if(i_spin == i_spin_add)
              {
                MLCluster::Spin const spin = {cluster.origin.site, -i_spin->pos};
                shifted.push_back(spin);
              }
              else
              {
                MLCluster::Spin const spin = {i_spin_add->site, i_spin_add->pos-i_spin->pos};
                shifted.push_back(spin);
              }
            t_uint const n_prime = cluster_to_integer(shifted, flavors, neighbors);
            if( n_prime < Nmax ) database[n_prime] = false;
            else is_in_shell = false;
          } // end of loop over spins.
        } // end of loop over cluster in newly added class.
        if( is_in_shell ) _out.push_back( clusters );
      } // end of loop over database.
    } // end of create_clusters.

  } // end of namespace CE
} // namespace LaDa

//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <algorithm>

#include <atat/vectmac.h>
#include <crystal/lattice.h>
#include <crystal/neighbors.h>
#include <crystal/symmetry_operator.h>
#include <crystal/which_site.h>
#include <opt/types.h>
#include <opt/fuzzy.h>

#include "create_pairs.h"


namespace LaDa
{
  namespace CE 
  {

    void create_pairs( const Crystal :: Lattice &_lat,
                       types::t_unsigned _max_neigh,
                       t_ClusterClasses &_out,
                       size_t _site )
    {
      __DEBUGTRYBEGIN
      typedef std::vector< std::vector<Cluster> > t_EquivClusters;
      if( _max_neigh == 0 ) return;

      // Creates a typical cluster.
      Cluster cluster;
      cluster.vectors.resize(2);
      cluster.vectors[0] = atat::rVector3d(0,0,0);
      atat::rVector3d const &origin(_lat.sites[_site].pos);
      cluster.eci = 0e0;

      size_t const N(_out.size()+_max_neigh);
      size_t size_try;
      switch( _max_neigh ) // includes enough neighbors for fcc or bcc. Hence for any lattice?
      {
        default: size_try = 30*_max_neigh; break;
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


      // Computes point group at lattice point.
      std::vector<Crystal::SymmetryOperator> point_group;
      foreach(Crystal::SymmetryOperator const &op, _lat.space_group)
        if( Crystal::which_site(op(origin), !_lat.cell, _lat.sites) == _site )
          point_group.push_back(op);
      std::vector<Crystal::SymmetryOperator>::const_iterator i_op;
      std::vector<Crystal::SymmetryOperator>::const_iterator const i_op_end = point_group.end();
      
      // goes over neighbors and adds to class.
      // Crystal::Neighbors gives us the insurance that classes which are found are complete,
      // since Crystal::Neighbors includes all neighbors up to a cutoff size.
      Crystal::Neighbors neighbors(size_try, origin);
      types::t_int first_of_size(_out.size()-1);
      Crystal::Neighbors::const_iterator i_first = neighbors.begin(_lat);
      Crystal::Neighbors::const_iterator const i_end = neighbors.end();
      types::t_real current_norm(0);
      for(; i_first != i_end; ++i_first)
      {
        // Checks if norm has changed.
        if( not Fuzzy::is_zero(i_first->distance - current_norm) )
        {
          current_norm = i_first->distance;
          ++first_of_size;
        }

        // checks if is in known class.
        atat::rVector3d pos(i_first->pos+origin);
        std::vector<Cluster> *clusters = NULL;
        bool not_found = true;
        for(size_t i(first_of_size); not_found and i < _out.size(); ++i)
        {
          clusters = &_out[i];
          t_Clusters::const_iterator i_begin;
          t_Clusters::const_iterator const i_end( clusters->end() );
          for(i_op = point_group.begin(); not_found and i_op != i_op_end; ++i_op)
          {
            atat::rVector3d const vec( (*i_op)(i_first->pos+origin) - (*i_op)(origin) );
            LADA_ASSERT( which_site( (*i_op)(pos), !_lat.cell, _lat.sites) != -1, "error" );
            LADA_ASSERT( which_site(vec + origin, !_lat.cell, _lat.sites) != -1, "error" );
            for(i_begin = clusters->begin(); not_found and i_begin != i_end; ++i_begin)
              if( Fuzzy::is_zero( atat::norm2(vec - i_begin->vectors[1]) ) ) not_found = false;
          }
        } 

        // creates new class if not found.
        if( not_found ) 
        {
          _out.resize(_out.size() + 1);
          clusters = &_out.back();
        }

        cluster.vectors[1] = i_first->pos;
        clusters->push_back(cluster);
      }
      if( _out.size() < N ) 
      {
        _out.clear();
        create_pairs( _lat, _max_neigh + 5, _out, _site);
      }
      _out.resize(N);

      __DEBUGTRYEND(,"Could not create pair clusters.\n" )
    }
   
  } // namespace CE
}

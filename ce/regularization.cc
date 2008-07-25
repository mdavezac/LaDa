//
//  Version: $Id: functional_builder.impl.h 685 2008-07-22 17:20:39Z davezac $
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "regularization.h"

namespace CE
{
  void find_pis( const std::vector< std::vector< Cluster > > &_clusters,
                 const Crystal::Structure & _str,
                 std::vector< types::t_int > &_pis )
  {
    typedef std::vector< Cluster >  t_Clusters;
    typedef std::vector< t_Clusters > t_EClusters;
    typedef std::vector< types::t_int > t_Pis;

    _pis.resize( _clusters.size() );

    atat::rMatrix3d inv_cell = !(~str.cell);
    Crystal :: Structure :: t_Atoms :: const_iterator i_atom = str.atoms.begin();
    Crystal :: Structure :: t_Atoms :: const_iterator i_atom_end = str.atoms.end();
    for(; i_atom != i_atom_end; ++i_atom) // loop over atoms
    {

      // loop over classes of clusters.
      t_EClusters :: const_iterator i_clusters = _clusters.begin();
      t_EClusters :: const_iterator i_clusters_end = _clusters.end();
      t_Pis :: iterator i_pi = pis.begin();
      for( ; i_clusters != i_clusters_end; ++i_clusters, ++i_pi ) 
      {
        *i_pi = 0;
        t_Clusters :: const_iterator i_cluster = i_clusters->begin();
        t_Clusters :: const_iterator i_cluster_end = i_clusters->end();
        // loop over equivalent clusters.
        for( ; i_cluster != i_cluster_end; ++i_cluster )
        {
          // null cluster
          if ( i_cluster->vectors.size() == 0 ) continue;

          // loop over cluster origin.
          typedef std::vector<atat::rVector3d> :: iterator vec_iterator;
          vec_iterator i_cpos_begin = i_cluster->vectors.begin();
          vec_iterator i_cpos_center = i_cluster->vectors.begin();
          vec_iterator i_cpos_end = i_cluster->vectors.end();
          vec_iterator i_cpos;
          for (; i_cpos_center != i_cpos_end; ++i_cpos_center ) 
          {   
            // loop over lattice positions in cluster.
            types::t_int result(1);
            for ( i_cpos = i_cpos_begin; i_cpos != i_cpos_end; ++i_cpos )
            {
              atat::rVector3d shift = i_atom->pos - *i_cpos_center;
              
              if (is_int( (!lattice->cell)*shift)) 
              {
                // finds atom to which lattice site is equivalent
                Crystal::Structure::t_Atoms::const_iterator
                  i_equiv = _str.atoms.begin()
                for (; i_equiv != i_atom_end; ++i_equiv)  
                  if ( atat::equivalent_mod_cell(*i_cpos + shift,
                                                 i_equiv->pos,inv_cell) ) 
                    break;
                
                __ASSERT( index == str.atoms.size(),
                            "Could not find equivalent of atom "
                          << (*i_cpos + shift)
                          << " in Lamarck::generate_functionals\n" )
                result *= i_equiv->type > 0 ? 1: -1;
              }

            }  // end of loop over cluster points
            *i_pi += result;
        
          } // end of rotation
        
        }  // end of loop over equivalent clusters
      } // end of loop over cluster classes.
    }  // end of loop over atoms 
  
    (*polynome) *= 1.0 /( (types::t_real) str.atoms.size() );
  
    // now computes constituent strain 
    return std::pair< t_Chemical*, t_CS* >( polynome, new t_CS(str) );
  }

  void find_pis( const std::vector< std::vector< Cluster > > &_clusters,
                 std::vector< Crystal::Structure > & _str;
                 std::vector< std::vector< types::t_int > > &_pis )
  {
    namespace bl = boost::lambda;
    _pis.resize( _str.size() );
    void (*ptr_func)( const std::vector< std::vector< Cluster > > &_clusters,
                      const Crystal::Structure & _str;
                      std::vector< types::t_int > &_pis ) = &find_pis;
    opt::concurrent_loop
    (
      _str.begin(), _str.end(), _pis.begin(),
      bl::bind( ptr_func, bl::constant( _clusters ), bl::_1, bl::_2 )
    );
  }
} // end of namespace CE


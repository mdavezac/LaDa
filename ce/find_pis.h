//
//  Version: $Id$
//
#ifndef LADA_CE_FIND_PIS_H
#define LADA_CE_FIND_PIS_H

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/lambda/bind.hpp>

#include<opt/algorithms.h>

#include "cluster.h"

namespace LaDa
{
  namespace CE
  {
    //! \brief Computes pis of \a _str for \a _clusters.
    //! \param[in] _cluster is a vector of containers of symmetrically equivalent
    //!                     cluster, centered on the origin.
    //! \param[in] _str the structure for which to compute the pis.
    //! \param[out] _pis the computed pis, one per class of symmetrically
    //!                  equivalent clusters.
    template< class T_CLUSTERS, class T_PIS >
      void find_pis( const T_CLUSTERS &_clusters,
                     const Crystal::Structure & _str,
                     T_PIS &_pis );
    //! \brief Computes pis of \a _str for \a _clusters.
    //! \see[in] _cluster is a vector of containers of symmetrically equivalent
    //!                     cluster, centered on the origin.
    //! \param[in] _str structures for which to compute the pis.
    //! \param[out] _pis the computed pis, one per class of symmetrically
    //!                  equivalent clusters.
    template< class T_CLUSTERS, class T_PIS >
      void find_pis( const T_CLUSTERS &_clusters,
                     const std::vector< Crystal::Structure > & _str,
                     T_PIS &_pis );
    template< class T_CLUSTERS, class T_PIS >
      void find_pis( const T_CLUSTERS &_clusters,
                     const Crystal::Structure & _str,
                     T_PIS &_pis )
      {
        namespace bl = boost::lambda;
        typedef T_CLUSTERS  t_Clusters;
        typedef typename t_Clusters :: value_type t_EClusters;
        typedef T_PIS t_Pis;
      
        _pis.resize( _clusters.size() );
        std::fill( _pis.begin(), _pis.end(), 0 );
      
        atat::rMatrix3d inv_cell = !(_str.cell);
        Crystal :: Structure :: t_Atoms :: const_iterator i_atom = _str.atoms.begin();
        Crystal :: Structure :: t_Atoms :: const_iterator i_atom_end = _str.atoms.end();
        for(; i_atom != i_atom_end; ++i_atom) // loop over atoms
        {
          // loop over classes of clusters.
          typename t_Clusters :: const_iterator i_clusters = _clusters.begin();
          typename t_Clusters :: const_iterator i_clusters_end = _clusters.end();
          typename t_Pis :: iterator i_pi = _pis.begin();
      
          for( ; i_clusters != i_clusters_end; ++i_clusters, ++i_pi ) 
          {
            types::t_real nm = 1e0 / types::t_real(   _str.atoms.size()
                                                    * i_clusters->size() );
            __ASSERT( i_clusters->size() == 0, "Cluster class is empty.\n" );
            if( not i_clusters->front().vectors.size() )
              { *i_pi = 1; continue; }
            typename t_EClusters :: const_iterator i_cluster = i_clusters->begin();
            typename t_EClusters :: const_iterator i_cluster_end = i_clusters->end();
            // loop over equivalent clusters.
            for( ; i_cluster != i_cluster_end; ++i_cluster )
            {
              // loop over cluster origin.
              typedef std::vector<atat::rVector3d> :: const_iterator vec_iterator;
              vec_iterator i_cpos_begin = i_cluster->vectors.begin();
              vec_iterator i_cpos_center = i_cluster->vectors.begin();
              vec_iterator i_cpos_end = i_cluster->vectors.end();
              vec_iterator i_cpos;
              for (; i_cpos_center != i_cpos_end; ++i_cpos_center ) 
              {   
                // loop over lattice positions in cluster.
                types::t_real result(1);
                for ( i_cpos = i_cpos_begin; i_cpos != i_cpos_end; ++i_cpos )
                {
                  atat::rVector3d shift = i_atom->pos - *i_cpos_center;
                  
                  if ( not is_int( (!_str.lattice->cell)*shift) ) continue;
                  
                  // finds atom to which lattice site is equivalent
                  Crystal::Structure::t_Atoms::const_iterator i_equiv = _str.atoms.begin();
                  for (; i_equiv != i_atom_end; ++i_equiv)  
                    if ( atat::equivalent_mod_cell( *i_cpos + shift, i_equiv->pos,inv_cell) ) 
                      break;
      
                  __ASSERT( i_equiv == i_atom_end,
                            "Could not find equivalent site.\n" )
                  result *= i_equiv->type;
      
                }  // end of loop over cluster points
                *i_pi += result * nm / (types::t_real) i_cluster->vectors.size();
              } // end of rotation
            }  // end of loop over equivalent clusters
          } // end of loop over cluster classes.
        }  // end of loop over atoms 
      }
      
    template< class T_CLUSTERS, class T_PIS >
      void find_pis( const T_CLUSTERS &_clusters,
                     const std::vector< Crystal::Structure > & _str,
                     T_PIS &_pis )
      {
        namespace bl = boost::lambda;
        typedef typename T_PIS::value_type t_Pis;
        _pis.resize( _str.size() );
        std::vector< Crystal::Structure > :: const_iterator i_str = _str.begin();
        std::vector< Crystal::Structure > :: const_iterator i_str_end = _str.end();
        typename T_PIS :: iterator i_pis = _pis.begin();
        for( ; i_str != i_str_end; ++i_str, ++i_pis )
          find_pis( _clusters, *i_str, *i_pis );
      }

  } // end of namespace CE
} // namespace LaDa
#endif

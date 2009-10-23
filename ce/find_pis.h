//
//  Version: $Id$
//
#ifndef LADA_CE_FIND_PIS_H
#define LADA_CE_FIND_PIS_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/numeric/ublas/vector.hpp>

#include <vector> 
#include <algorithm> 

#include <crystal/lattice.h>
#include <crystal/structure.h>
#include <crystal/which_site.h>
#include <crystal/smith.h>


#include "cluster.h"

namespace LaDa
{
  namespace CE
  {
    //! \brief Computes the pis of a structure.
    //! \param[in]  _cls classes of equivalent clusters
    //! \param[in]  _str structure for which to compute the pis.
    //! \param[out] _out Pis on output.
    //! \param[in]  _site sublattice for which to compute the pis.
    //!                   If negative or larger than the number of sublattices,
    //!                   clusters will be computed for all sublattices.
    template< class T_CLUSTERCLASSES, class T_VECTOR >
      void find_pis( T_CLUSTERCLASSES const &_cls,
                     Crystal::Structure const &_str, 
                     T_VECTOR &_out, 
                     types::t_int _site = -1) 
      {
        LADA_ASSERT( _str.lattice != NULL, "lattice not set.\n" )
        LADA_ASSERT( _str.lattice->sites.size() != 0, "lattice does not contain sites.\n" )
      

        bool const single_site( _site >= 0 );
        size_t const site( std::abs(_site) );
        _out.resize(_cls.size());
        std::fill(_out.begin(), _out.end(), 0);
        if( _cls.size() == 0 ) return;
      
        std::vector< std::vector<size_t> > atomic_map;
        Crystal::get_smith_map( _str, atomic_map );
        bool const dosite( _str.lattice->sites.size() != 1 );
        size_t const Npersite( _str.atoms.size() / _str.lattice->sites.size() );
        atat::rMatrix3d const inv_cell( !_str.lattice->cell );
        Crystal::Lattice::t_Sites const &sites(_str.lattice->sites);
      
        Crystal::t_SmithTransform const
          transform( Crystal::get_smith_transform(_str.lattice->cell, _str.cell) );
      
        Crystal::Structure::t_Atoms::const_iterator i_first = _str.atoms.begin();
        Crystal::Structure::t_Atoms::const_iterator const i_end = _str.atoms.end();
        for(; i_first != i_end; ++i_first) // loop over atomic positions.
        {
          if( single_site and i_first->site != site ) continue;
      
          typename T_CLUSTERCLASSES :: const_iterator i_class = _cls.begin();
          typename T_CLUSTERCLASSES :: const_iterator i_class_end = _cls.end();
          typename T_VECTOR :: iterator i_out = _out.begin();
          for(; i_class != i_class_end; ++i_class, ++i_out) // loop over cluster classes.
          {
            size_t const Nperclass( i_class->size() );
            if( Nperclass == 0 ) continue; // nothin in class.

            typedef typename T_CLUSTERCLASSES :: value_type :: const_iterator t_cit;
            t_cit i_cluster = i_class->begin();
            size_t const order( i_cluster->size() );
            if( order == 0 ) { *i_out = 1e0; continue; } // zero order class.

            t_cit const i_cluster_end = i_class->end();
            types::t_real pi(0);
            for(; i_cluster != i_cluster_end; ++i_cluster ) // loop over equivalent clusters.
            {
              typedef std::vector<atat::rVector3d> :: const_iterator t_cit;
              t_cit i_vec_begin = i_cluster->vectors.begin();
              t_cit const i_vec_end = i_cluster->vectors.end();
              t_cit i_center = i_vec_begin;
              for(; i_center != i_vec_end; ++i_center ) // loop over cluster centers.
              {
                types::t_real fig(1);
                std::vector<atat::rVector3d> :: const_iterator i_vec = i_vec_begin;
                atat::rVector3d const shift( i_first->pos - *i_center);
                if( Crystal::which_site(shift, inv_cell, sites) == -1 ) continue;
                for(; i_vec != i_vec_end; ++i_vec) // loop over cluster spins.
                {
                  atat::rVector3d const site_pos(*i_vec + shift);
                  types::t_int const sindex
                  ( 
                    dosite ? Crystal::which_site(site_pos, inv_cell, sites): 0
                  );
                  LADA_DOASSERT( sindex != -1, "Site not found.\n" )
                  LADA_DOASSERT( sindex < _str.lattice->sites.size(), "Index out of range.\n" );
                  size_t const site_index(sindex);
                  atat::rVector3d const pos( site_pos - sites[site_index].pos );
                  size_t const smith_index = Crystal::get_linear_smith_index(transform, pos);
                  fig *= _str.atoms[ atomic_map[site_index][smith_index] ].type;
                } // end of loop over spins.
        
                pi += fig;
              }
            } // loop over equivalent clusters.
            *i_out += pi / types::t_real(Npersite*Nperclass*order);
          } // loop over classes of clusters.
        } // loop over atomic positions.
#       ifdef LADA_DEBUG
          foreach(types::t_real &o, _out) if( Fuzzy::is_zero(o) ) o = 0e0;
#       endif
      }

    //! \brief Computes pis of \a _str for \a _clusters.
    //! \param[in] _cluster is a vector of containers of symmetrically equivalent
    //!                     cluster, centered on the origin.
    //! \param[in] _str the structure for which to compute the pis.
    //! \param[out] _pis the computed pis, one per class of symmetrically
    //!                  equivalent clusters.
    template< class T_CLUSTERS, class T_PIS >
      void find_pis( const T_CLUSTERS &_clusters,
                     const std::vector< Crystal::Structure > & _str,
                     T_PIS &_pis, size_t _site = 0 )
      {
        typedef typename T_PIS::value_type t_Pis;
        _pis.resize( _str.size() );
        std::vector< Crystal::Structure > :: const_iterator i_str = _str.begin();
        std::vector< Crystal::Structure > :: const_iterator i_str_end = _str.end();
        typename T_PIS :: iterator i_pis = _pis.begin();
        for( ; i_str != i_str_end; ++i_str, ++i_pis )
          find_pis( _clusters, *i_str, *i_pis, _site );
      }

    template< class T_CLUSTERS, class T_PIS >
    void find_pis2( const T_CLUSTERS &_clusters,
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
  } // end of namespace CE
} // namespace LaDa
#endif

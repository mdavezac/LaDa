//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <crystal/symmetry_operator.h>
#include <opt/types.h>
#include <opt/fuzzy.h>

#include "find_pis.h"


namespace LaDa
{
  namespace CE 
  {
    template< class T_VECTOR >
    void find_pis_impl( t_ClusterClasses const &_cls,
                        Crystal::Structure const &_str, 
                        T_VECTOR &_out, 
                        size_t _site )
    {
      LADA_ASSERT( _str.lattice != NULL, "lattice not set.\n" )

      _out.resize(_cls.size());
      std::fill(_out.begin(), _out.end(), 0);
      if( _cls.size() == 0 ) return;

      std::vector< std::vector<size_t> > atomic_map;
      Crystal::get_smith_map( _str, atomic_map );
      bool const dosite( _str.lattice->sites.size() != 1 );
      atat::rMatrix3d const inv_cell( !_str.lattice->cell );
      Crystal::Lattice::t_Sites const &sites(_str.lattice->sites);

      Crystal::t_SmithTransform const transform( Crystal::get_smith_transform(_str.lattice->cell, _str.cell) );

      Crystal::Structure::t_Atoms::const_iterator i_first = _str.atoms.begin();
      Crystal::Structure::t_Atoms::const_iterator const i_end = _str.atoms.end();
      for(; i_first != i_end; ++i_first)
      {
        if( i_first->site != _site ) continue;

        t_ClusterClasses :: const_iterator i_class = _cls.begin();
        t_ClusterClasses :: const_iterator const i_class_end = _cls.end();
        typename T_VECTOR :: iterator i_out = _out.begin();
        for(; i_class != i_class_end; ++i_class, ++i_out)
        {
          if( i_class->size() == 0 ) continue;
          t_Clusters :: const_iterator i_cluster = i_class->begin();
          t_Clusters :: const_iterator const i_cluster_end = i_class->end();
          types::t_real fig(1);
          for(; i_cluster != i_cluster_end; ++i_cluster ) 
          {
            std::vector<atat::rVector3d> :: const_iterator i_vec = i_cluster->vectors.begin();
            std::vector<atat::rVector3d> :: const_iterator const i_vec_end = i_cluster->vectors.end();
            for(; i_vec != i_vec_end; ++i_vec)
            {
              atat::rVector3d const site_pos(*i_vec + i_first->pos);
              types::t_int const sindex
              ( 
                dosite ? Crystal::which_site(site_pos, inv_cell, sites): 0
              );
              LADA_DOASSERT( sindex != -1, "Site not found.\n" );
              size_t const site_index(sindex);
              atat::rVector3d const pos( site_pos - sites[site_index].pos );
              size_t const smith_index = Crystal::get_linear_smith_index(transform, pos);
              fig *= _str.atoms[ atomic_map[site_index][smith_index] ].type;
            }

            *i_out += fig;
          };
        }
      }
    }



    void find_pis2( t_ClusterClasses const &_cls,
                    const Crystal::Structure &_str, 
                    std::vector<types::t_real> &_out, 
                    size_t _site)
      { find_pis_impl( _cls, _str, _out, _site ); }

    void find_pis2( t_ClusterClasses const &_cls,
                    const Crystal::Structure &_str, 
                    boost::numeric::ublas::vector<types::t_real> &_out, 
                    size_t _site )
      { find_pis_impl( _cls, _str, _out, _site ); }
  } // namespace CE
}

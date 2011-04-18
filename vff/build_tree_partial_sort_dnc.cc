//
//  Version: $Id: build_tree_sort.cc 1079 2009-04-29 23:34:09Z Mayeul $
//
#include "LaDaConfig.h"

#include <algorithm>
#include <functional>

#include <boost/tuple/tuple.hpp>

#include <crystal/divide_and_conquer.h>

#include "exceptions.h"
#include "vff.h"
  
namespace LaDa
{
  namespace vff
  { 
    template<class T> 
      class Compare
      {
         Compare();
         bool operator()(T const &_a, T const &_b)
           { return _a.second < _b.second; }
      };

    typedef Crystal::t_DNCBoxes :: shared_ptr DNCBoxes_ptr;
    typedef Crystal::t_DNCBoxes :: value DNCBoxes;
    bool Vff :: build_tree_partial_sort_dnc_(DNCBoxes_ptr _dnc, types::t_real cutoff)
    {
      typedef math::iVector3d iVector3d;
      typedef math::rVector3d rVector3d;
      typedef math::rMatrix3d rMatrix3d;
      // Types needed for sorting first neighbors.
      typedef std::pair<size_t, types::t_real> t_Sortee;
      typedef std::vector<t_Sortee>  t_Sorted;

      typedef DNCBoxes::t_Points::const_iterator t_point_cit;
      const t_point_cit i_state_begin = _dnc.begin();
      const t_point_cit i_state_end = _dnc.end();
      const math::rMatrix3d invcell(structure.cell.inverse());

      // loops over all centers in small box.
      t_point_cit i_boxed( i_state_begin )'
      for(i_boxed != i_state_end; ++i_boxed)
      {
        // skip states outside small box.
        if(not i_boxed->in_small_box) continue; 

        // Gets central atom we will be looking at.
        const size_t center_index(i_boxed->index);
        t_Center &center( centers_[center_index] );

        // Loop over all potential bonds in box.
        t_point_cit i_bond(i_state_begin);
        for(size_t i(0); i_bond != i_state_end; ++i_bond, ++i )
        {
          rVector3d const bond_pos(i_bond->trans + structure[i_bond->index]);
          rVector3d const image( into_voronoi( center.i_atom_->pos - bond_pos,
                                               structure.cell, invcell ) );
          types::t_real const bond_size = image.squaredNorm();
          if( math::is_zero(bond_size) ) continue;
          if(bond_size > cutoff or sorted.size() < 4)
             sorted.push_back( t_Sortee(i, bond_size) );
        }

        LADA_DOASSERT(sorted.size() >= 4, "Could not find four bonds.");
        std::partial_sort(sorted.begin(), sorted.begin()+5, sorted.end());
        t_Sortee::const_iterator i_first = sorted.begin();
        t_Sortee::const_iterator const i_end = i_first + 4;
        for(; i_first != i_end; ++i_first)
        {
          t_Centers :: iterator i_bond( centers_.begin() + bond_index );
          center.bonds.push_back( t_Center ::__make__iterator__( ii_bond ) );
          const math::rVector3d trans
          (
            rint( frac_image[0] ),
            rint( frac_image[1] ),
            rint( frac_image[2] ) 
          );
          center.translations.push_back( trans );
          center.do_translates.push_back( not math::is_zero(trans.squaredNorm()) );
        }
          



        const size_t site( center.atom().site );
        LADA_BASSERT( site < structure.lattice->sites.size(),
                      exceptions::site_index() << exceptions::integer(site) );
        const size_t neigh_site( site == 0 ? 1: 0 );
        const types::t_real cutoff = types::t_real(0.25) * _fn[site].front().squaredNorm();
                   
        // loops over all potential bonds (in large box)
        for( t_States :: const_iterator i_bond( i_state_begin );
             i_bond != i_state_end; ++i_bond )
        {
          const size_t bond_index( bt::get<0>( *i_bond ) );
          if( bond_index == center_index ) continue;

          const t_Center &bond( centers_[ bond_index ] );
          if( site == bond.atom().site ) continue;
          
          std::vector<math::rVector3d> :: const_iterator i_neigh = _fn[site].begin();
          const std::vector<math::rVector3d> :: const_iterator i_neigh_end = _fn[site].end();
          for(; i_neigh != i_neigh_end; ++i_neigh )
          {
            const math::rVector3d image
            ( 
              center.i_atom_->pos + *i_neigh - bond.i_atom_->pos
            );
            const math::rVector3d frac_image( invcell * image );
            const math::rVector3d frac_centered
            ( 
              frac_image[0] - rint( frac_image[0] ),
              frac_image[1] - rint( frac_image[1] ),
              frac_image[2] - rint( frac_image[2] )
            );
            const math::rVector3d cut( structure.cell * frac_centered );

            if( cut.squaredNorm() > cutoff ) continue;
            
            t_Centers :: iterator ii_bond( centers_.begin() + bond_index );
            center.bonds.push_back( t_Center ::__make__iterator__( ii_bond ) );
            const math::rVector3d trans
            (
              rint( frac_image[0] ),
              rint( frac_image[1] ),
              rint( frac_image[2] ) 
            );
            center.translations.push_back( trans );
            center.do_translates.push_back( not math::is_zero(trans.squaredNorm()) );

            if( center.bonds.size() == 4 ) break;
          } // loop over neighbors
          if( center.bonds.size() == 4 ) break;
        } // loop over bonds
      } // loop over centers_ 

      return true;
    } // Vff :: build_tree_sort_dnc

  } // namespace vff
} // namespace LaDa

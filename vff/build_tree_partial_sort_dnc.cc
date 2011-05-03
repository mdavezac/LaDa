#include "LaDaConfig.h"

#include <algorithm>
#include <functional>

#include <boost/tuple/tuple.hpp>

#include <crystal/periodic_dnc.h>
#include <crystal/lattice.h>

#include "exceptions.h"
#include "vff.h"
  
namespace LaDa
{
  namespace vff
  { 
    template<class T> 
      struct Compare
      {
         Compare() {};
         bool operator()(T const &_a, T const &_b) const
           { return _a.second < _b.second; }
      };

    bool Vff :: build_tree_partial_sort_dnc( Crystal::DnCBoxes::value_type const & _dnc,
                                             types::t_real cutoff )
    {
      typedef math::iVector3d iVector3d;
      typedef math::rVector3d rVector3d;
      typedef math::rMatrix3d rMatrix3d;
      // Types needed for sorting first neighbors.
      typedef std::pair<size_t, types::t_real> t_Sortee;
      typedef std::vector<t_Sortee>  t_Sorted;

      typedef Crystal::DnCBoxes::value_type::const_iterator t_point_cit;
      const t_point_cit i_state_begin = _dnc.begin();
      const t_point_cit i_state_end   = _dnc.end();
      const math::rMatrix3d invcell(structure.cell.inverse());

      // loops over all centers in small box.
      t_point_cit i_boxed( i_state_begin );
      for(; i_boxed != i_state_end; ++i_boxed)
      {
        // skip states outside small box.
        if(not i_boxed->in_small_box) continue; 

        // Position of current atom being investigated.
        rVector3d current_pos = structure.atoms[i_boxed->index].pos + i_boxed->translation;

        // Gets central atom we will be looking at.
        const size_t center_index(i_boxed->index);
        t_Center &center( centers_[center_index] );

        // Loop over all potential bonds in box.
        t_point_cit i_bond(i_state_begin);
        t_Sorted sorted;
        for(size_t i(0); i_bond != i_state_end; ++i_bond, ++i )
        {
          rVector3d const bond_pos(structure.atoms[i_bond->index].pos + i_bond->translation);
          types::t_real const bond_size = (bond_pos - current_pos).squaredNorm(); 
          if( bond_size > types::tolerance and (bond_size < cutoff or sorted.size() < 4) )
             sorted.push_back( t_Sortee(i, bond_size) );
        }

        LADA_DOASSERT(sorted.size() >= 4, "Could not find four bonds.");
        std::partial_sort(sorted.begin(), sorted.begin()+4, sorted.end(), Compare<t_Sortee>());
        t_Sorted::const_iterator i_first = sorted.begin();
        t_Sorted::const_iterator const i_end = i_first + 4;
        for(; i_first != i_end; ++i_first)
        {
          Crystal::DnCBoxes::value_type::value_type const &point = _dnc[i_first->first];
          t_Centers :: iterator ii_bond( centers_.begin() + point.index );
          center.bonds.push_back( t_Center ::__make__iterator__( ii_bond ) );
          const math::rVector3d
            trans( math::round(invcell * (point.translation - i_boxed->translation)));
          center.translations.push_back( trans );
          center.do_translates.push_back( not math::is_zero(trans.squaredNorm()) );
        }

      } // loop over points in box.

      return true;
    } // Vff :: build_tree_sort_dnc

  } // namespace vff
} // namespace LaDa

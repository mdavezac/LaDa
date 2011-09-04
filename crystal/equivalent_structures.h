#ifndef LADA_CRYSTAL_EQUIVALENT_STRUCTURES_H
#define LADA_CRYSTAL_EQUIVALENT_STRUCTURES_H
#include "LaDaConfig.h"

#include <list>
#include <vector>
#include <algorithm>

#include <boost/foreach.hpp>

#include <math/misc.h>
#include <math/fuzzy.h>
#include <math/gruber.h>

#include "structure.h"
#include "utilities.h"
#include "space_group.h"


namespace LaDa
{
  namespace crystal 
  {
    namespace details
    {
      template<class T_TYPE> 
        int min_test(TemplateStructure<T_TYPE> const &_a)
        {
          std::vector<size_t> mini(_a.size(), 1);
          std::vector<size_t>::iterator i_fmin = mini.begin();
          typename TemplateStructure<T_TYPE>::const_iterator i_first = _a.begin();
          typename TemplateStructure<T_TYPE>::const_iterator const i_end = _a.end();
          for(; i_first != i_end; ++i_first, ++i_fmin)
          {
            if(*i_fmin > 1) continue;
            CompareOccupations<T_TYPE> cmp(i_first->type);
            typename TemplateStructure<T_TYPE>::const_iterator i_second = i_first + 1;
            std::vector<size_t>::iterator i_smin = i_fmin + 1;
            for(; i_second != i_end; ++i_second, ++i_smin)
              if(cmp(i_second->type)) ++(*i_fmin), ++(*i_smin);
          }
          return std::min_element(mini.begin(), mini.end()) - mini.begin();
        }
    }
    //! \brief Returns true if two structures are equivalent. 
    //! \details Two structures are equivalent in a crystallographic sense,
    //!          e.g. without reference to cartesian coordinates or possible
    //!          motif rotations which leave the lattice itself invariant. A
    //!          supercell is *not* equivalent to its lattice, unless it is a
    //!          trivial supercell.
    //! \param[in] _a: The first structure.
    //! \param[in] _b: The second structure.
    //! \param[in] with_scale: whether to take the scale into account. Defaults to true.
    //! \param[in] tolerance: Tolerance when comparing distances. Defaults to
    //!            types::t_real. It is in the same units as the structures scales, if
    //!            that is taken into account, otherwise, it is in the same
    //!            units as _a.scale.
    template<class T_TYPE> 
      bool equivalent( TemplateStructure<T_TYPE> const &_a,
                       TemplateStructure<T_TYPE> const &_b,
                       bool with_scale = true,
                       types::t_real _tol = types::tolerance )
      {
        // different number of atoms.
        if(_a.size() != _b.size()) return false;
        types::t_real const scaleA = _a.scale();
        types::t_real const scaleB
          = with_scale ?
              _b.scale():
              _a.scale() * std::pow(_a.cell().determinant() / _b.cell().determinant(), 1e0/3e0);
        // different volume.
        if(math::neq( (_a.cell()*scaleA).determinant(),
                      (_b.cell()*scaleB).determinant(), 3e0*_tol)) return false;
        
        // check possible rotation. 
        math::rMatrix3d const cellA = math::gruber(_a.cell(), 100, _tol) * scaleA;
        math::rMatrix3d const cellB = math::gruber(_b.cell(), 100, _tol) * scaleB;
        math::rMatrix3d const invA = cellA.inverse();
        math::rMatrix3d const invB = cellB.inverse();
        math::rMatrix3d const rot = cellA * cellB.inverse();
        if(not math::is_identity(rot * (~rot), 2*_tol)) return false;
        if(math::neq(rot.determinant(), 1e0, 3*_tol))  return false;
        
        // Now checks atomic sites. 
        // first computes point-group symmetries.
        boost::shared_ptr<t_SpaceGroup> pg = cell_invariants(cellA);
        // then find the occupation type with the smallest number of occurences.
        typename TemplateStructure<T_TYPE>::const_reference minatom = _a[details::min_test(_a)];
        CompareOccupations<T_TYPE> mincheck(minatom.type);
        
        // Computes possible translations, looking at only one type of site-occupation.
        // The center of gravity will tell us about a possible translation of
        // the cartesian basis.
        math::rVector3d transA(0,0,0);
        size_t nA(0);
        foreach(typename TemplateStructure<T_TYPE>::const_reference atomA, _a)
          if(mincheck(atomA.type))
          {
            transA += into_voronoi(atomA.pos * scaleA, cellA, invA);
            ++nA;
          }
        transA /= types::t_real(nA);
        math::rVector3d transB(0,0,0);
        size_t nB = 0;
        foreach(typename TemplateStructure<T_TYPE>::const_reference atomB, _b)
          if(mincheck(atomB.type))
          {
            transB += into_voronoi(atomB.pos * scaleB, cellB, invB);
            ++nB;
            if(nB > nA) return false;
          }
        transB /= types::t_real(nB);

        // loop over possible motif rotations.
        foreach(math::Affine3d const &invariant, *pg)
        {
          // creates a vector referencing B atomic sites.
          // Items from this list will be removed as they are found.
          typedef std::list<size_t> t_List;
          t_List atomsA;
          for(size_t i(0); i < _a.size(); ++i) atomsA.push_back(i);

          math::rMatrix3d const rotation = invariant.linear() * rot;
          
          typename TemplateStructure<T_TYPE>::const_iterator i_b = _b.begin();
          typename TemplateStructure<T_TYPE>::const_iterator const i_bend = _b.end();
          for(; i_b != i_bend; ++i_b)
          {
            math::rVector3d const pos = rotation * (into_voronoi(i_b->pos*scaleB, cellB, invB) - transB);
            CompareOccupations<T_TYPE> const cmp(i_b->type);
            typename t_List :: iterator i_first =  atomsA.begin();
            typename t_List :: iterator const i_end = atomsA.end();
            if(i_first == i_end) return false;
            for(; i_first != i_end; ++i_first)
            {
              if( not math::is_integer(invA * (pos - _a[*i_first].pos*scaleA + transA), 4*_tol) ) continue;
              if( cmp(_a[*i_first].type) ) break;
            }
            if(i_first == i_end) break;
            atomsA.erase(i_first);
          }
          if(i_b == i_bend) return true;
        }
        return false;
      }
    
  } // namespace crystal
} // namespace LaDa
#endif

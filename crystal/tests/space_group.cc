#include "LaDaConfig.h"
#include <iostream>
#include <string>
#include <vector>

#include <boost/foreach.hpp>

#include <opt/debug.h>
#include "../space_group.h"
#include "../structure.h"
#include "../equivalent_structures.h"
#include "../primitive.h"

namespace LaDa
{
  namespace crystal
  {
    //! \brief Returns true if two structures are equivalent lattices in same cartesian coordinates. 
    //! \details Two structures are equivalent in a mathematical sense, for the
    //!          same cartesian coordinates. That is, all sites defined by one
    //!          lattice correspond to the same site in the other lattice.
    //!          Supercells are *not* equivalent to their parent lattice (use
    //!          in conjunction with primitive() to get this behavior). Note that since the carete
    //! \param[in] _a: The first structure.
    //! \param[in] _b: The second structure.
    //! \param[in] with_scale: whether to take the scale into account. Defaults to true.
    //! \param[in] tolerance: Tolerance when comparing distances. Defaults to
    //!            types::t_real. It is in the same units as the structures scales, if
    //!            that is taken into account, otherwise, it is in the same
    //!            units as _a.scale.
    bool equivalent_lattices( Structure< LADA_TYPE > const &_a,
                              Structure< LADA_TYPE > const &_b,
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
      // different lattice parameterization.
      if(not math::is_integer(_a.cell().inverse() * _b.cell(), 3e0*_tol)) return false;
      if(not math::is_integer(_b.cell().inverse() * _a.cell(), 3e0*_tol)) return false;
      
      // check possible rotation. 
      math::rMatrix3d const cellA = math::gruber(_a.cell(), 100, _tol) * scaleA;
      math::rMatrix3d const invA = cellA.inverse();
      
      // creates a vector referencing A atomic sites.
      // Items from this list will be removed as they are found.
      typedef std::list<size_t> t_List;
      t_List atomsA;
      for(size_t i(0); i < _a.size(); ++i) atomsA.push_back(i);

      Structure< LADA_TYPE >::const_iterator i_b = _b.begin();
      Structure< LADA_TYPE >::const_iterator const i_bend = _b.end();
      for(; i_b != i_bend; ++i_b)
      {
        math::rVector3d const pos = into_voronoi(i_b->pos()*scaleB, cellA, invA);
        CompareOccupations< LADA_TYPE > const cmp(i_b->type());
        t_List :: iterator i_first =  atomsA.begin();
        t_List :: iterator const i_end = atomsA.end();
        if(i_first == i_end) return false;
        for(; i_first != i_end; ++i_first)
        {
          if( not math::is_integer(invA * (pos - _a[*i_first]->pos*scaleA), 4*_tol) ) continue;
          if( cmp(_a[*i_first]->type) ) break;
        }
        if(i_first == i_end) break;
        atomsA.erase(i_first);
      }
      if(i_b == i_bend) return true;
      return false;
    }

    Structure< LADA_TYPE > b5(types::t_real u)
    {
      Structure< LADA_TYPE > lattice;
      types::t_real const x(u), y(0.25 - u);
      lattice.set_cell(0, 0.5, 0.5)
                      (0.5, 0, 0.5)
                      (0.5, 0.5, 0);
      lattice.add_atom(5.000000e-01, 5.000000e-01, 5.000000e-01, "A")
                      (5.000000e-01, 2.500000e-01, 2.500000e-01, "A")
                      (2.500000e-01, 5.000000e-01, 2.500000e-01, "A")
                      (2.500000e-01, 2.500000e-01, 5.000000e-01, "A")
                      (8.750000e-01, 8.750000e-01, 8.750000e-01, "B")
                      (1.250000e-01, 1.250000e-01, 1.250000e-01, "B")
                      (     x,     x,     x, "X")
                      (     x,     y,     y, "X")
                      (     y,     x,     y, "X")
                      (     y,     y,     x, "X")
                      (    -x,    -x,    -x, "X")
                      (    -x,    -y,    -y, "X")
                      (    -y,    -x,    -y, "X")
                      (    -y,    -y,    -x, "X");
      return lattice;
    }
  }
}


int main()
{
  using namespace LaDa;
  using namespace LaDa::crystal;
  using namespace LaDa::math;
  Structure< LADA_TYPE > lattice;
  lattice.set_cell(0, 0.5, 0.5)
                  (0.5, 0, 0.5)
                  (0.5, 0.5, 0);
  lattice.add_atom(0,0,0, "Si");
  boost::shared_ptr<t_SpaceGroup> sg = space_group(lattice);
  LADA_DOASSERT(sg->size() == 48, "Improper number of rotations.");
  foreach(t_SpaceGroup::const_reference aff, *sg)
  {
    LADA_DOASSERT(math::is_null(aff.translation()), "Non-zero translation.\n");
    Structure<LADA_TYPE> transformed = lattice.transform(aff);
    equivalent_lattices(lattice, transformed);
  }
  
  lattice = b5(0.25);
  sg = space_group(lattice);
  LADA_DOASSERT(sg->size() == 48, "Improper number of rotations.");
  foreach(t_SpaceGroup::const_reference aff, *sg)
  {
    Structure<LADA_TYPE> transformed = lattice.transform(aff);
    equivalent_lattices(lattice, transformed);
  }

  lattice = b5(0.36);
  sg = space_group(lattice);
  LADA_DOASSERT(sg->size() == 48, "Improper number of rotations.");
  foreach(t_SpaceGroup::const_reference aff, *sg)
  {
    Structure<LADA_TYPE> transformed = lattice.transform(aff);
    equivalent_lattices(lattice, transformed);
  }
  return 0;
}

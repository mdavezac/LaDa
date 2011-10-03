#include "LaDaConfig.h"
#include<iostream>
#include<string>
#include<vector>

#include <boost/foreach.hpp>
#include "../space_group.h"
#include <opt/debug.h>

int main()
{
  using namespace LaDa;
  using namespace LaDa::crystal;
  using namespace LaDa::math;
  rMatrix3d cell;
  cell << 0, 0.5, 0.5, 0.5, 0, 0.5, 0.5, 0.5, 0;
  boost::shared_ptr<t_SpaceGroup> invariants = cell_invariants(cell);
  LADA_DOASSERT(invariants->size() == 48, "Improper number of rotations.");
  foreach(math::Affine3d const &aff, *invariants)
  {
    LADA_DOASSERT(math::is_null(aff.translation()), "Non-zero translation.\n");
    math::rMatrix3d newcell = aff * cell;
    LADA_DOASSERT(math::is_integer(cell * newcell.inverse()), "Not invariant.\n");
    LADA_DOASSERT(math::is_integer(newcell * cell.inverse()), "Not invariant.\n");
  }

  cell << -0.5,0.5,0.5, 0.5,-0.5,0.5, 0.5,0.5,-0.5;
  invariants = cell_invariants(cell);
  LADA_DOASSERT(invariants->size() == 48, "Improper number of rotations.");
  foreach(math::Affine3d const &aff, *invariants)
  {
    LADA_DOASSERT(math::is_null(aff.translation()), "Non-zero translation.\n");
    math::rMatrix3d newcell = aff * cell;
    LADA_DOASSERT(math::is_integer(cell * newcell.inverse()), "Not invariant.\n");
    LADA_DOASSERT(math::is_integer(newcell * cell.inverse()), "Not invariant.\n");
  }

  cell << -0.6,0.5,0.5, 0.6,-0.5,0.5, 0.6,0.5,-0.5;
  invariants = cell_invariants(cell);
  LADA_DOASSERT(invariants->size() == 4, "Improper number of rotations.");
  foreach(math::Affine3d const &aff, *invariants)
  {
    LADA_DOASSERT(math::is_null(aff.translation()), "Non-zero translation.\n");
    math::rMatrix3d newcell = aff * cell;
    LADA_DOASSERT(math::is_integer(cell * newcell.inverse()), "Not invariant.\n");
    LADA_DOASSERT(math::is_integer(newcell * cell.inverse()), "Not invariant.\n");
  }

  cell << -0.7,0.7,0.7, 0.6,-0.5,0.5, 0.6,0.5,-0.5;
  invariants = cell_invariants(cell);
  LADA_DOASSERT(invariants->size() == 8, "Improper number of rotations.");
  foreach(math::Affine3d const &aff, *invariants)
  {
    LADA_DOASSERT(math::is_null(aff.translation()), "Non-zero translation.\n");
    math::rMatrix3d newcell = aff * cell;
    LADA_DOASSERT(math::is_integer(cell * newcell.inverse()), "Not invariant.\n");
    LADA_DOASSERT(math::is_integer(newcell * cell.inverse()), "Not invariant.\n");
  }

  cell << -0.765,0.7,0.7, 0.665,-0.5,0.5, 0.6,0.5,-0.5;
  invariants = cell_invariants(cell);
  LADA_DOASSERT(invariants->size() == 2, "Improper number of rotations.");
  foreach(math::Affine3d const &aff, *invariants)
  {
    LADA_DOASSERT(math::is_null(aff.translation()), "Non-zero translation.\n");
    math::rMatrix3d newcell = aff * cell;
    LADA_DOASSERT(math::is_integer(cell * newcell.inverse()), "Not invariant.\n");
    LADA_DOASSERT(math::is_integer(newcell * cell.inverse()), "Not invariant.\n");
  }
  return 0;
}

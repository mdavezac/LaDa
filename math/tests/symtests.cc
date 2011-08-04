#include "LaDaConfig.h"

#include<iostream>

#include <opt/debug.h>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/serialization.hpp>

#include "../serialize.h"
#include "../misc.h"
#include "../symmetry_operator.h"

using namespace std;
int main()
{
  using namespace LaDa;
  using namespace LaDa::math;

  Affine3d eig;
  eig = AngleAxis(0, rVector3d::UnitX());
  LADA_ASSERT(is_isometry(eig), "Is not an isometry.\n");
  LADA_ASSERT(is_rotation(eig), "Is not a rotation.\n");
  eig = AngleAxis(0, rVector3d::UnitX()) * Translation(0,1,0);
  LADA_ASSERT(not is_rotation(eig), "Is a rotation.\n");
  LADA_ASSERT(is_isometry(eig), "Is not an isometry.\n");
  eig.matrix()(0,0) = 5;
  LADA_ASSERT(not is_isometry(eig), "Is an isometry.\n");
  LADA_ASSERT(not is_rotation(eig), "Is not a rotation.\n");
  eig.matrix()(0,1) = 1;
  eig.linear() /= eig.linear().determinant();
  LADA_ASSERT(not is_isometry(eig), "Is an isometry.\n");
  eig = AngleAxis(2.*pi/3., rVector3d(1,1,1));
  LADA_ASSERT(not is_isometry(eig), "Is an isometry.\n");
  LADA_ASSERT(not is_rotation(eig), "Is not a rotation.\n");

  rMatrix3d mat;
  mat << -0.5, 0.5, 0.5, 0.5, -0.5, 0.5, 0.5, 0.5, -0.5;
  eig = AngleAxis(2.*pi/3., rVector3d(1,1,1) / std::sqrt(3.)) * Translation(0,1.5,0); 
  LADA_ASSERT(is_invariant(eig, mat), "Is not invariant.\n");
  eig = AngleAxis(pi/3., rVector3d(1,1,1) / std::sqrt(3.)) * Translation(0,1.5,0); 
  LADA_ASSERT(not is_invariant(eig, mat), "Is invariant.\n");

  rVector3d vec(1,1,1);
  eig = AngleAxis(pi/2., rVector3d::UnitZ());
  LADA_ASSERT( not is_invariant(eig, vec), "Is invariant");
  eig = Translation(2, 0, 0) * eig;
  LADA_ASSERT( is_invariant(eig, vec), "Is not invariant");
  return 0;
}

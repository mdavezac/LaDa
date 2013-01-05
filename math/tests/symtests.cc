#include "PyladaConfig.h"

#include<iostream>

#define PYLADA_DOASSERT(a,b) \
        { \
          if((not (a)))\
          { \
            std::cerr << __FILE__ << ", line: " << __LINE__ << "\n" << b; \
            throw 0;\
          }\
        }

#include "../math.h"

using namespace std;
int main()
{
  using namespace Pylada;
  using namespace Pylada::math;

  Affine3d eig;
  eig = AngleAxis(0, rVector3d::UnitX());
  PYLADA_DOASSERT(is_isometry(eig), "Is not an isometry.\n");
  PYLADA_DOASSERT(is_rotation(eig), "Is not a rotation.\n");
  eig = AngleAxis(0, rVector3d::UnitX()) * Translation(0,1,0);
  PYLADA_DOASSERT(not is_rotation(eig), "Is a rotation.\n");
  PYLADA_DOASSERT(is_isometry(eig), "Is not an isometry.\n");
  eig.matrix()(0,0) = 5;
  PYLADA_DOASSERT(not is_isometry(eig), "Is an isometry.\n");
  PYLADA_DOASSERT(not is_rotation(eig), "Is not a rotation.\n");
  eig.matrix()(0,1) = 1;
  eig.linear() /= eig.linear().determinant();
  PYLADA_DOASSERT(not is_isometry(eig), "Is an isometry.\n");
  eig = AngleAxis(2.*pi/3., rVector3d(1,1,1));
  PYLADA_DOASSERT(not is_isometry(eig), "Is an isometry.\n");
  PYLADA_DOASSERT(not is_rotation(eig), "Is not a rotation.\n");

  rMatrix3d mat;
  mat << -0.5, 0.5, 0.5, 0.5, -0.5, 0.5, 0.5, 0.5, -0.5;
  eig = AngleAxis(2.*pi/3., rVector3d(1,1,1) / std::sqrt(3.)) * Translation(0,1.5,0); 
  PYLADA_DOASSERT(is_invariant(eig, mat), "Is not invariant.\n");
  eig = AngleAxis(pi/3., rVector3d(1,1,1) / std::sqrt(3.)) * Translation(0,1.5,0); 
  PYLADA_DOASSERT(not is_invariant(eig, mat), "Is invariant.\n");

  rVector3d vec(1,1,1);
  eig = AngleAxis(pi/2., rVector3d::UnitZ());
  PYLADA_DOASSERT( not is_invariant(eig, vec), "Is invariant");
  eig = Translation(2, 0, 0) * eig;
  PYLADA_DOASSERT( is_invariant(eig, vec), "Is not invariant");
  return 0;
}

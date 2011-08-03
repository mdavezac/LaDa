#include "LaDaConfig.h"

#include<iostream>

#include <opt/debug.h>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/serialization.hpp>

#include "../serialize.h"
#include "../misc.h"

using namespace std;
int main()
{
  using namespace LaDa;
  using namespace LaDa::math;

  Affine3d eig, eig2, eig3;
  eig =   AngleAxis( 0.3,  rVector3d::UnitZ())
        * AngleAxis(-0.15, rVector3d::UnitX())
        * AngleAxis( 1.15, rVector3d::UnitZ())
        * Translation(0, 15, 0.5);
  eig =   AngleAxis(0, rVector3d::UnitZ())
        * AngleAxis(0, rVector3d::UnitX())
        * AngleAxis(0, rVector3d::UnitZ())
        * Translation(0, 0, 0);
  std::ostringstream ss;
  boost::archive::text_oarchive oa( ss );
  oa << eig;

  std::istringstream ss2( ss.str().c_str() );
  boost::archive::text_iarchive ia( ss2 );
  ia >> eig2;

  LADA_DOASSERT(math::eq(eig, eig2), "Could not reload eig.\n");
  
  rVector3d vec(1, 0, 0);
  eig =   AngleAxis( pi / 2,  rVector3d::UnitZ())
        * Translation(0, 0, 0);
  LADA_ASSERT( eq(eig*vec, rVector3d(0, 1, 0)), "Did not rotate as anticipated.\n");
  eig *= AngleAxis( pi / 2,  rVector3d::UnitY());
  LADA_ASSERT( eq(eig*vec, rVector3d(0, 0, -1)), "Did not rotate as anticipated.\n");
  eig = Translation(0, 1.5, 0) * eig; 
  eig2 = eig;
  eig3 = eig;
  LADA_ASSERT( eq(eig*vec, rVector3d(0, 1.5, -1)), "Did not rotate as anticipated.\n");
  eig *= AngleAxis( -pi / 2,  rVector3d::UnitZ());
  LADA_ASSERT( eq(eig*vec, rVector3d(1, 1.5,  0)), "Did not rotate as anticipated.\n");
  eig2 *= AngleAxis( -pi / 2,  rVector3d::UnitX());
  LADA_ASSERT( eq(eig2*vec, rVector3d(0, 1.5,  -1)), "Did not rotate as anticipated.\n");
  eig3 *= AngleAxis( -pi / 2,  rVector3d::UnitX());
  LADA_ASSERT( eq(eig3*vec, rVector3d(0, 1.5,  -1)), "Did not rotate as anticipated.\n");

  return 0;
}

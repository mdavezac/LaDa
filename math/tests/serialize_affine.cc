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
  eig =   AngleAxis( pi / 2,  rVector3d::UnitX())
        * Translation(0, 0, 0);
  LADA_ASSERT( eq(eig*vec, rVector3d(1, 0, 0)), "Did not rotate as anticipated.\n");
  eig = AngleAxis( pi / 2,  rVector3d::UnitY());
  LADA_ASSERT( eq(eig*vec, rVector3d(0, 0, -1)), "Did not rotate as anticipated.\n");
  eig = AngleAxis( pi / 2,  rVector3d::UnitZ());
  LADA_ASSERT( eq(eig*vec, rVector3d(0, 1, 0)), "Did not rotate as anticipated.\n");
  eig = AngleAxis( pi / 2,  rVector3d::UnitZ()) * Translation(0,0.5,0);
  LADA_ASSERT( eq(eig*vec, rVector3d(-0.5, 1, 0)), "Did not rotate as anticipated.\n");
  eig = AngleAxis( pi / 2,  rVector3d::UnitZ());
  eig = Translation(0,0.5,0) * eig;
  LADA_ASSERT( eq(eig*vec, rVector3d(0, 1.5, 0)), "Did not rotate as anticipated.\n");
  LADA_ASSERT( eq(eig, eig), "Did not compare as expected."); 
  LADA_ASSERT( neq(eig, eig * Translation(0,-0.5,0)), "Did not compare as expected.");
  LADA_ASSERT( neq(eig, eig * AngleAxis(pi/4, rVector3d::UnitX())), "Did not compare as expected.");
  LADA_ASSERT( eq(eig, eig * AngleAxis(types::tolerance, rVector3d::UnitX()), 8. * types::tolerance),
               "Did not compare as expected.");
  LADA_ASSERT( neq(eig, eig * AngleAxis(types::tolerance, rVector3d::UnitX()), types::tolerance),
               "Did not compare as expected.");
  return 0;
}

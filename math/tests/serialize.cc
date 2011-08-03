#include "LaDaConfig.h"

#include<iostream>

#include <opt/debug.h>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/serialization.hpp>

#include "../serialize.h"
#include "../misc.h"

#if LADA_TEST_INCTYPE == 0
#  define LADA_TYPE  Eigen::Matrix<types::t_real, LADA_X, LADA_Y>
#elif LADA_TEST_INCTYPE == 1
#  define LADA_TYPE  Eigen::Matrix<types::t_int, LADA_X, LADA_Y>
#elif LADA_TEST_INCTYPE == 6
#  define LADA_TYPE  LaDa::math::rVector3d
#  define LADA_X 3
#  define LADA_Y 1
#elif LADA_TEST_INCTYPE == 7
#  define LADA_TYPE  LaDa::math::iVector3d
#  define LADA_X 3
#  define LADA_Y 1
#elif LADA_TEST_INCTYPE == 8
#  define LADA_TYPE  LaDa::math::rMatrix3d
#  define LADA_X 3
#  define LADA_Y 3
#elif LADA_TEST_INCTYPE == 9
#  define LADA_TYPE  LaDa::math::iMatrix3d
#  define LADA_X 3
#  define LADA_Y 3
#endif

using namespace std;
int main()
{
  using namespace LaDa;
  using namespace LaDa::math;

  LADA_TYPE eig, eig2;
  for(size_t i(0); i < LADA_X; ++i)
    for(size_t j(0); j < LADA_Y; ++j)
    {
      eig(i, j) = i + j * LADA_X;
      eig2(i, j) = 0;
    }
  std::ostringstream ss;
  boost::archive::text_oarchive oa( ss );
  oa << eig;

  std::istringstream ss2( ss.str().c_str() );
  boost::archive::text_iarchive ia( ss2 );
  ia >> eig2;

  LADA_DOASSERT(math::eq(eig, eig2), "Could not reload eig.\n");
  return 0;
}

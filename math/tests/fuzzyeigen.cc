#include "LaDaConfig.h"

#include<iostream>

#define LADA_DOASSERT(a,b) \
        { \
          if((not (a)))\
          { \
            std::cerr << __FILE__ << ", line: " << __LINE__ << "\n" << b; \
            throw 0;\
          }\
        }

#include "../misc.h"

#if LADA_TEST_INCTYPE == 0
#  define LADA_UNITY math::iVector3d(1,1,1)
#  define LADA_ZERO math::iVector3d(0,0,0)
#  define LADA_SMALL math::iVector3d(0, 0, 0)
#  define LADA_TOL   1
#elif LADA_TEST_INCTYPE == 1
#  define LADA_ZERO math::iMatrix3d::Zero()
#  define LADA_SMALL math::iMatrix3d::Zero()
#  define LADA_TOL   1
#elif LADA_TEST_INCTYPE == 2
#  define LADA_UNITY math::rVector3d(1,1,1)
#  define LADA_ZERO math::rVector3d(0,0,0)
#  define LADA_SMALL math::rVector3d(0.5 * types::tolerance, 0, 0)
#  define LADA_TOL   types::tolerance
#elif LADA_TEST_INCTYPE == 3
#  define LADA_ZERO math::rMatrix3d::Zero()
#  define LADA_TOL   types::tolerance
#endif

using namespace std;
int main()
{
  using namespace LaDa;
  using namespace LaDa::math;
# if LADA_TEST_INCTYPE == 1
    iMatrix3d LADA_UNITY;
    LADA_UNITY << 1, 0, 0, 0, 1, 0, 0, 0, 1;
# elif LADA_TEST_INCTYPE == 3
    rMatrix3d LADA_UNITY, LADA_SMALL;
    LADA_UNITY << 1, 0, 0, 0, 1, 0, 0, 0, 1;
    LADA_SMALL << 0.5 * types::tolerance, 0, 0, 0, 0.5 * types::tolerance, 0, 0.5 * types::tolerance, 0, 0;
# endif

  LADA_DOASSERT(is_identity(LADA_UNITY+LADA_SMALL), "unexpected unity.")
# if LADA_TEST_TYPE > 1
  LADA_DOASSERT(not is_identity(LADA_UNITY+10*LADA_SMALL, LADA_TOL), "unexpected unity.")
# endif
  LADA_DOASSERT(not is_identity(2*LADA_UNITY+LADA_SMALL, LADA_TOL), "unexpected unity.")
  LADA_DOASSERT(is_null(LADA_ZERO + LADA_SMALL), "unexpected null.")
  LADA_DOASSERT(not is_null(LADA_UNITY + LADA_SMALL, LADA_TOL), "unexpected null.")
  LADA_DOASSERT(is_integer(LADA_UNITY + LADA_SMALL), "unexpected integer.")
  LADA_DOASSERT(eq(LADA_UNITY + LADA_SMALL, LADA_UNITY), "unexpected eq.")
  LADA_DOASSERT(neq(2*LADA_UNITY, LADA_UNITY), "unexpected neq.")
# if LADA_TEST_TYPE > 1
    LADA_DOASSERT(not is_integer(LADA_UNITY + LADA_SMALL, LADA_TOL), "unexpected integer.")
# endif
# if LADA_TEST_TYPE == 3
    matrix << 0, 0.5, 0.5, 0.5, 0, 0.5, 0.5, 0.5, 0;
    LADA_DOASSERT(are_periodic_images( rVector3d(0.5, 0.25, 0.5),
                                     rVector3d(0.5, 0.25, 0.5) + matrix * rVector3d(1, -10, 2),
                                     matrix.inverse() ), "NOT PERIODIC.\n" );
# endif
  return 0;
}

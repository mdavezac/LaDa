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

#include "../fuzzy.h"

#if LADA_TEST_INCTYPE == 0
#  define LADA_FIRST 5
#  define LADA_SMALL 0
#  define LADA_TOL   1
#elif LADA_TEST_INCTYPE == 1
#  define LADA_FIRST 5e0
#  define LADA_SMALL   0.5 * types::tolerance
#  define LADA_TOL   10 * types::tolerance
#endif

using namespace std;
int main()
{
  using namespace LaDa;
  using namespace LaDa::math;

  LADA_DOASSERT(eq(LADA_FIRST, LADA_FIRST+LADA_SMALL), "unexpected eq result.\n")
  LADA_DOASSERT(not eq(LADA_FIRST, LADA_FIRST+1), "unexpected eq result.\n")
  LADA_DOASSERT(not lt(LADA_FIRST, LADA_FIRST+LADA_SMALL), "unexpected lt result.\n")
  LADA_DOASSERT(not lt(LADA_FIRST, LADA_FIRST+LADA_SMALL, LADA_TOL), "unexpected lt result.\n")
  LADA_DOASSERT(lt(LADA_FIRST, LADA_FIRST+5), "unexpected lt result.\n")
  LADA_DOASSERT(not gt(LADA_FIRST, LADA_FIRST-LADA_SMALL), "unexpected gt result.\n")
  LADA_DOASSERT(not gt(LADA_FIRST, LADA_FIRST-LADA_SMALL, LADA_TOL), "unexpected gt result.\n")
  LADA_DOASSERT(gt(LADA_FIRST+5, LADA_FIRST), "unexpected gt result.\n")
  LADA_DOASSERT(is_identity(1 + LADA_SMALL), "unexpected unity.")
  LADA_DOASSERT(not is_identity(LADA_FIRST), "unexpected unity.")
  LADA_DOASSERT(is_null(LADA_SMALL), "unexpected null.")
  LADA_DOASSERT(not is_null(LADA_FIRST), "unexpected null.")
  LADA_DOASSERT(is_integer(LADA_FIRST + LADA_SMALL), "unexpected integer.")
# if LADA_TEST_TYPE != 0
    LADA_DOASSERT(not is_integer(LADA_FIRST + LADA_SMALL, LADA_TOL), "unexpected integer.")
# endif
  return 0;
}

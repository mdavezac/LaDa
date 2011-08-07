#include "LaDaConfig.h"
#include<iostream>
#include<string>
#include<vector>

#include "../compare_sites.h"
#if LADA_INC_TYPE == 0
#  define LADA_TYPE std::string
#  define LADA_INIT_TYPE0 a.type = "Au";
#  define LADA_INIT_TYPE1 b.type = "Pd";
#elif LADA_INC_TYPE == 1
#  define LADA_TYPE int
#  define LADA_INIT_TYPE0 a.type = 0;
#  define LADA_INIT_TYPE1 b.type = 2;
#elif LADA_INC_TYPE == 2
#  define LADA_TYPE types::t_real
#  define LADA_INIT_TYPE0 a.type = 0.0;
#  define LADA_INIT_TYPE1 b.type = 2.0;
#elif LADA_INC_TYPE == 3
#  define LADA_TYPE std::vector<std::string>
#  define LADA_INIT_TYPE0 a.type.push_back("Au"); a.type.push_back("Pd");
#  define LADA_INIT_TYPE1 b.type.push_back("Au"); b.type.push_back("Pd");
#endif

using namespace std;
int main()
{
  using namespace LaDa;
  using namespace LaDa::crystal;
  using namespace LaDa::math;
  Atom< LADA_TYPE > a, b;
  a.pos = rVector3d(0.1, 0.2, 0.1);
  b.pos = rVector3d(0.3, 0.2, 0.1);
  LADA_INIT_TYPE0;
  LADA_INIT_TYPE1;
  LADA_ASSERT(not compare_sites(a)(b.type), "equivalent.\n");
  LADA_ASSERT(not compare_sites(a)(b.pos), "equivalent.\n");
  LADA_ASSERT(not compare_sites(a)(b), "equivalent.\n");
  LADA_ASSERT(compare_sites(a)(a.pos), "not equivalent.\n");
  LADA_ASSERT(compare_sites(a)(a.type), "not equivalent.\n");
  LADA_ASSERT(compare_sites(a)(a), "not equivalent.\n");
  b.pos = a.pos;
  LADA_ASSERT(not compare_sites(a)(b), "not equivalent.\n");
  b.pos += rVector3d::UnitX();
  b.type = a.type;
  LADA_ASSERT(not compare_sites(a)(b), "not equivalent.\n");
  b.pos = a.pos;
  LADA_ASSERT(compare_sites(a)(b), "equivalent.\n");

  return 0;
}

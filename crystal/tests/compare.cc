#include "LaDaConfig.h"
#include<iostream>
#include<string>
#include<vector>

#include "../compare_sites.h"
#if LADA_TEST_INCTYPE == 0
#  define LADA_TYPE std::string
#  define LADA_INIT_TYPE0 a.type = "Au";
#  define LADA_INIT_TYPE1 b.type = "Pd";
#elif LADA_TEST_INCTYPE == 1
#  define LADA_TYPE int
#  define LADA_INIT_TYPE0 a.type = 0;
#  define LADA_INIT_TYPE1 b.type = 2;
#elif LADA_TEST_INCTYPE == 2
#  define LADA_TYPE types::t_real
#  define LADA_INIT_TYPE0 a.type = 0.0;
#  define LADA_INIT_TYPE1 b.type = 2.0;
#elif LADA_TEST_INCTYPE == 3
#  define LADA_TYPE std::vector<std::string>
#  define LADA_INIT_TYPE0 a.type.push_back("Au"); a.type.push_back("Pd");
#  define LADA_INIT_TYPE1 b.type.push_back("Au");
#elif LADA_TEST_INCTYPE == 4
#  define LADA_TYPE std::vector<int>
#  define LADA_INIT_TYPE0 a.type.push_back(0); a.type.push_back(5);
#  define LADA_INIT_TYPE1 b.type.push_back(0);
#elif LADA_TEST_INCTYPE == 5
#  define LADA_TYPE std::vector<types::t_real>
#  define LADA_INIT_TYPE0 a.type.push_back(0e0); a.type.push_back(5e0);
#  define LADA_INIT_TYPE1 b.type.push_back(0e0);
#elif LADA_TEST_INCTYPE == 6
#  define LADA_TYPE std::set<std::string>
#  define LADA_INIT_TYPE0 a.type.insert("Au"); a.type.insert("Pd");
#  define LADA_INIT_TYPE1 b.type.insert("Au");
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
  std::cout << a << "\n";
  LADA_DOASSERT(not compare_sites(a)(b.type), "equivalent.\n");
  LADA_DOASSERT(not compare_sites(a)(b.pos), "equivalent.\n");
  LADA_DOASSERT(not compare_sites(a)(b), "equivalent.\n");
  LADA_DOASSERT(compare_sites(a)(a.pos), "not equivalent.\n");
  LADA_DOASSERT(compare_sites(a)(a.type), "not equivalent.\n");
  LADA_DOASSERT(compare_sites(a)(a), "not equivalent.\n");
  b.pos = a.pos;
  LADA_DOASSERT(not compare_sites(a)(b), "equivalent.\n");
  LADA_DOASSERT(compare_sites(a)(b.pos), "not equivalent.\n");
  b.pos += rVector3d::UnitX();
  b.type = a.type;
  LADA_DOASSERT(not compare_sites(a)(b), "equivalent.\n");
  LADA_DOASSERT(compare_sites(a)(b.type), "not equivalent.\n");
  b.pos = a.pos;
  LADA_DOASSERT(compare_sites(a)(b), "equivalent.\n");
# if LADA_TEST_INCTYPE >= 3 and LADA_TEST_INCTYPE <= 5
  b.type.pop_back();
  LADA_DOASSERT(not compare_sites(a)(b), "equivalent.\n");
  LADA_DOASSERT(not compare_sites(b)(a), "equivalent.\n");
# endif 

  return 0;
}

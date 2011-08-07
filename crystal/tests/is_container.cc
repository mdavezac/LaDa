
#include "LaDaConfig.h"
#include<iostream>
#include<string>
#include<vector>

#include "../compare_sites.h"
#if LADA_TEST_INCTYPE == 0
#  define LADA_TYPE std::string
#  define LADA_TEST false
#  define LADA_INIT_TYPE var = "Au";
#  define LADA_STRING_TEST "Au"
#elif LADA_TEST_INCTYPE == 1
#  define LADA_TYPE int
#  define LADA_TEST false
#  define LADA_INIT_TYPE var = 0;
#  define LADA_STRING_TEST "0"
#elif LADA_TEST_INCTYPE == 2
#  define LADA_TYPE types::t_real
#  define LADA_TEST false
#  define LADA_INIT_TYPE var = 0.0;
#  define LADA_STRING_TEST "0"
#elif LADA_TEST_INCTYPE == 3
#  define LADA_TYPE std::vector<std::string>
#  define LADA_TEST true
#  define LADA_INIT_TYPE var.push_back("Au"); var.push_back("Pd");
#  define LADA_STRING_TEST "Au, Pd"
#elif LADA_TEST_INCTYPE == 4
#  define LADA_TYPE std::vector<int>
#  define LADA_TEST true
#  define LADA_INIT_TYPE var.push_back(0); var.push_back(5);
#  define LADA_STRING_TEST "0, 5"
#elif LADA_TEST_INCTYPE == 5
#  define LADA_TYPE std::vector<types::t_real>
#  define LADA_TEST true
#  define LADA_INIT_TYPE var.push_back(0e0); var.push_back(5e0);
#  define LADA_STRING_TEST "0, 5"
#endif

using namespace std;
int main()
{
  using namespace LaDa;
  using namespace LaDa::crystal;
  using namespace LaDa::math;
  LADA_ASSERT(details::is_container<LADA_TYPE>::value == LADA_TEST, "Did not complete test.\n");
  LADA_TYPE var;
  LADA_INIT_TYPE;
  LADA_ASSERT(details::print_occupation(var) == LADA_STRING_TEST, "Did not print.\n");

  return 0;
}



#include "LaDaConfig.h"
#include<iostream>
#include<string>
#include<vector>
#include "../structure.h"
#include "../is_container.h"
#include <opt/debug.h>

#if LADA_TEST_INCTYPE == 0
# define TYPE std::string
# define TYPE0 "Au"
# define TYPE1 "Pd"
# define TEST_TYPE \
    LADA_DOASSERT(structure[0].type == "Au", "structure atom data incorrect.\n"); \
    LADA_DOASSERT(structure[1].type == "Pd", "structure atom data incorrect.\n");
#elif LADA_TEST_INCTYPE == 1
# define TYPE std::vector< std::string >
# define TYPE0 "Au"
# define TYPE1 "Au", "Pd"
# define TEST_TYPE \
    LADA_DOASSERT(structure[0].type.size() == 1, "Incorrect number of types in site 0."); \
    LADA_DOASSERT(structure[1].type.size() == 2, "Incorrect number of types in site 0."); \
    LADA_DOASSERT(structure[0].type[0] == "Au", "Incorrect type in site 0."); \
    LADA_DOASSERT(structure[1].type[0] == "Au", "Incorrect type in site 1."); \
    LADA_DOASSERT(structure[1].type[1] == "Pd", "Incorrect type in site 1."); 
#elif LADA_TEST_INCTYPE == 2
# define TYPE std::set< std::string >
  template<class T> TYPE create_set(T const &_t)
    { TYPE result; result.insert(_t); return result; }
  template<class T> TYPE create_set(T const &_t, T const &_t1)
    { TYPE result; result.insert(_t); result.insert(_t1); return result; }
# define TYPE0 "Au"
# define TYPE1 "Au", "Pd"
# define TEST_TYPE \
    LADA_DOASSERT(structure[0].type.size() == 1, "Incorrect number of types in site 0."); \
    LADA_DOASSERT(structure[1].type.size() == 2, "Incorrect number of types in site 0."); \
    LADA_DOASSERT(structure[0].type == create_set("Au"), "Incorrect type in site 0."); \
    LADA_DOASSERT(structure[1].type == create_set("Au", "Pd"), "Incorrect type in site 1."); 
#elif LADA_TEST_INCTYPE == 3
# define TYPE int
# define TYPE0 0
# define TYPE1 1
# define TEST_TYPE \
    LADA_DOASSERT(structure[0].type == 0, "structure atom data incorrect.\n"); \
    LADA_DOASSERT(structure[1].type == 1, "structure atom data incorrect.\n");
#elif LADA_TEST_INCTYPE == 4
# define TYPE std::vector<int>
# define TYPE0 0
# define TYPE1 0, 1
# define TEST_TYPE \
    LADA_DOASSERT(structure[0].type.size() == 1, "Incorrect number of types in site 0."); \
    LADA_DOASSERT(structure[1].type.size() == 2, "Incorrect number of types in site 0."); \
    LADA_DOASSERT(structure[0].type[0] == 0, "Incorrect type in site 0."); \
    LADA_DOASSERT(structure[1].type[0] == 0, "Incorrect type in site 1."); \
    LADA_DOASSERT(structure[1].type[1] == 1, "Incorrect type in site 1."); 
#elif LADA_TEST_INCTYPE == 5
# define TYPE std::set<int>
  template<class T> TYPE create_set(T const &_t)
    { TYPE result; result.insert(_t); return result; }
  template<class T> TYPE create_set(T const &_t, T const &_t1)
    { TYPE result; result.insert(_t); result.insert(_t1); return result; }
# define TYPE0 0
# define TYPE1 0, 1
# define TEST_TYPE \
    LADA_DOASSERT(structure[0].type.size() == 1, "Incorrect number of types in site 0."); \
    LADA_DOASSERT(structure[1].type.size() == 2, "Incorrect number of types in site 0."); \
    LADA_DOASSERT(structure[0].type == create_set(0), "Incorrect type in site 0."); \
    LADA_DOASSERT(structure[1].type == create_set(0, 1), "Incorrect type in site 1."); 
#elif LADA_TEST_INCTYPE == 6
# define TYPE LaDa::types::t_real
# define TYPE0 0
# define TYPE1 3.1416
# define TEST_TYPE \
    LADA_DOASSERT(structure[0].type == 0, "structure atom data incorrect.\n"); \
    LADA_DOASSERT(structure[1].type == 3.1416, "structure atom data incorrect.\n");
#elif LADA_TEST_INCTYPE == 7
# define TYPE std::vector<LaDa::types::t_real>
# define TYPE0 0
# define TYPE1 0, 3.1416
# define TEST_TYPE \
    LADA_DOASSERT(structure[0].type.size() == 1, "Incorrect number of types in site 0."); \
    LADA_DOASSERT(structure[1].type.size() == 2, "Incorrect number of types in site 0."); \
    LADA_DOASSERT(structure[0].type[0] == 0, "Incorrect type in site 0."); \
    LADA_DOASSERT(structure[1].type[0] == 0, "Incorrect type in site 1."); \
    LADA_DOASSERT(structure[1].type[1] == 3.1416, "Incorrect type in site 1."); 
#elif LADA_TEST_INCTYPE == 8
# define TYPE std::set<LaDa::types::t_real>
  template<class T> TYPE create_set(T const &_t)
    { TYPE result; result.insert(_t); return result; }
  template<class T> TYPE create_set(T const &_t, T const &_t1)
    { TYPE result; result.insert(_t); result.insert(_t1); return result; }
# define TYPE0 0
# define TYPE1 0, 3.1416
# define TEST_TYPE \
    LADA_DOASSERT(structure[0].type.size() == 1, "Incorrect number of types in site 0."); \
    LADA_DOASSERT(structure[1].type.size() == 2, "Incorrect number of types in site 0."); \
    LADA_DOASSERT(structure[0].type == create_set(0e0), "Incorrect type in site 0."); \
    LADA_DOASSERT(structure[1].type == create_set(0e0, 3.1416), "Incorrect type in site 1."); 
#endif


using namespace std;
int main()
{
  using namespace LaDa;
  using namespace LaDa::crystal;
  TemplateStructure< TYPE > structure;
  structure.set_cell(-0.5,0.5,0.5)
                    (0.5,-0.5,0.5)
                    (0.5,0.5,-0.5);
  structure.add_atom(0,0,0, TYPE0) 
                    (math::rVector3d::Ones() * 0.25, TYPE1);
  for(size_t i(0); i < 3; ++i)
  {
    for(size_t j(0); j < 3; ++j)
      if(i == j) { LADA_DOASSERT(structure(i,j) == -0.5, "structure cell data incorrect."); }
      else { LADA_DOASSERT(structure(i,j) == 0.5, "structure cell data incorrect."); }
    LADA_DOASSERT(structure[0].pos[i] == 0, "structure atom data incorrect.");
    LADA_DOASSERT(structure[1].pos[i] == 0.25, "structure atom data incorrect.");
  }
  TEST_TYPE;
  return 0;
}

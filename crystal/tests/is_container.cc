
#include "LaDaConfig.h"
#include <iostream>
#include <string>
#include <vector>

#include <opt/types.h>
#include <opt/debug.h>

#include "../is_container.h"

using namespace std;
template<class T> void check(T const &_var, bool _cont, bool _set, bool _string, 
                             std::string const &_str)
{
  using namespace LaDa::crystal::details;
  LADA_DOASSERT(is_container<T>::value == _cont, "");
  LADA_DOASSERT(is_set<T>::value == _set, "");
  LADA_DOASSERT(is_string<T>::value == _string, "");
  LADA_DOASSERT(is_nonstring_scalar<T>::value == not (_cont or _set or _string), "");
  LADA_DOASSERT(is_scalar<T>::value == not (_cont or _set), "");
  LADA_DOASSERT(print_occupation(_var) == _str, "");
}
int main()
{
  using namespace LaDa;
  using namespace LaDa::crystal;

  check<std::string>("Au", false, false, true, "Au");
  check<int>(1, false, false, false, "1");
  check<types::t_real>(3.1416, false, false, false, "3.1416");

  { std::vector<std::string> var(2); var[0] = "Au"; var[1] = "Pd"; 
    check(var, true, false, false, "Au, Pd"); }
  { std::vector<int> var(2); var[0] = 0; var[1] = 1; 
    check(var, true, false, false, "0, 1"); }
  { std::vector<types::t_real> var(2); var[0] = 0e0; var[1] = 3.1416; 
    check(var, true, false, false, "0, 3.1416"); }

  { std::set<std::string> var; var.insert("Au"); var.insert("Pd"); 
    check(var, false, true, false, "Au, Pd"); }
  { std::set<int> var; var.insert(0); var.insert(1); 
    check(var, false, true, false, "0, 1"); }
  { std::set<types::t_real> var; var.insert(0e0); var.insert(3.1416); 
    check(var, false, true, false, "0, 3.1416"); }

  return 0;
}

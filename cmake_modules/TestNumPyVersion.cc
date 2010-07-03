#include <iostream>
#include <Python.h>
//#include <numpy/__multiarray_api.h>
#include <numpy/ndarrayobject.h>

using namespace std;

int main ()
{
  int status (0);

#ifdef NPY_VERSION
  uint version = NPY_VERSION;
  std::cout << version << std::endl;
  return version;
#else
  return status;
#endif

}

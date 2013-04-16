#include <Python.h>
#include <numpy/ndarraytypes.h>

int main() {
   
  npy_intp dims[1] = {1};
  PyArrayObject *dummy = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
  PyArray_ENABLEFLAGS(dummy, 0);
    
  return 0;
}

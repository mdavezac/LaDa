#ifndef PYLADA_INCLUDE_NUMPY
// Makes sure that NPY_ARRAY style stuff exists, and ENABLEFLAGS
#if NUMPY_VERSION_MINOR >= 7
#  define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#endif
#include <numpy/arrayobject.h>
#ifndef LADA_NPY_NEWDEFS
#  define NPY_ARRAY_C_CONTIGUOUS NPY_C_CONTIGUOUS 
#  define NPY_ARRAY_WRITEABLE    NPY_WRITEABLE
#endif
#ifndef LADA_NPY_HAS_ENABLEFLAGS
#  define PyArray_ENABLEFLAGS(array, flags) array->flags |= flags
#  define PyArray_CLEARFLAGS(array, flags)  array->flags &= (!flags)
#  define PyArray_SetBaseObject(array, base) array->base = base
#endif

#endif

#include <Python.h>
#include <numpy/ndarrayobject.h>
#include <boost/mpl/int.hpp>

template<class T> class type;
      
template<> struct type<npy_longdouble> : public boost::mpl::int_<NPY_LONGDOUBLE> 
{
  typedef npy_longdouble np_type;
};
template<> struct type<npy_double> : public boost::mpl::int_<NPY_DOUBLE> 
{
  typedef npy_double np_type;
};
int main() {return 0;}

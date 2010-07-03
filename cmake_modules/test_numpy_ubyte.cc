#include <Python.h>
#include <numpy/ndarrayobject.h>
#include <numpy/npy_common.h>
#include <boost/mpl/int.hpp>


template<class T> class type;
template<> struct type<npy_ubyte> : public boost::mpl::int_<NPY_UBYTE> 
{
  typedef npy_ubyte np_type;
};
template<> struct type<npy_ubyte> : public boost::mpl::int_<NPY_BOOL> 
{
  typedef npy_bool np_type;
};
int main() {return type<npy_bool>::value;}

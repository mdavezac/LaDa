//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/python.hpp>

#include <opt/types.h>

#include "fxvector.h"
#include "vectmac.h"

#ifdef _TYPE_
#error "Please change _TYPE_ to something else, as it is already being used."
#endif
#ifdef _PYTHONNAME_
#error "Please change _PYTHONNAME_ to something else, as it is already being used."
#endif


namespace atat
{
  template<class T_TYPE>
  Vector3d<T_TYPE>* to_atatvector( boost::python::object & _ob )
  {
    T_TYPE x = boost::python::extract<T_TYPE>( _ob[0] );
    T_TYPE y = boost::python::extract<T_TYPE>( _ob[1] );
    T_TYPE z = boost::python::extract<T_TYPE>( _ob[2] );
    return new Vector3d<T_TYPE>(x,y,z);
  }
  template<class T_TYPE>
  Matrix3d<T_TYPE>* to_atatmatrix( boost::python::object & _ob )
  {
    boost::python::object a =  _ob[0];
    T_TYPE xx = boost::python::extract<T_TYPE>( a[0] );
    T_TYPE xy = boost::python::extract<T_TYPE>( a[1] );
    T_TYPE xz = boost::python::extract<T_TYPE>( a[2] );
    a =  _ob[1];
    T_TYPE yx = boost::python::extract<T_TYPE>( a[0] );
    T_TYPE yy = boost::python::extract<T_TYPE>( a[1] );
    T_TYPE yz = boost::python::extract<T_TYPE>( a[2] );
    a =  _ob[2];
    T_TYPE zx = boost::python::extract<T_TYPE>( a[0] );
    T_TYPE zy = boost::python::extract<T_TYPE>( a[1] );
    T_TYPE zz = boost::python::extract<T_TYPE>( a[2] );
    
    Matrix3d<T_TYPE> *m = new Matrix3d<T_TYPE>;
    m->x[0][0] = xx; m->x[0][1] = xy; m->x[0][2] = xz;
    m->x[1][0] = yx; m->x[1][1] = yy; m->x[1][2] = yz;
    m->x[2][0] = zx; m->x[2][1] = zy; m->x[2][2] = zz;
    return m;
  }
}



BOOST_PYTHON_MODULE(atat)
{
  using namespace boost::python;
  using namespace atat;
#  define _TYPE_ types::t_int
#  define _PYTHONNAME_(object) "i"#object
#  include "atat.py.hpp"

// #  define _TYPE_ types::t_real
// #  define _PYTHONNAME_(object) "r"+ #object
// #  include "atat.py.hpp"
}

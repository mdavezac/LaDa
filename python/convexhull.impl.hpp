//
//  Version: $Id$
//
#ifndef _PYTHONLADA_CONVEXHUL_IMPL_H_
#define _PYTHONLADA_CONVEXHUL_IMPL_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/python/object.hpp>
#include <string>

#include <opt/convex_hull.h>
#include <opt/types.h>
#include <opt/debug.h>

namespace PythonLaDa
{
  template< class T_TYPE > 
    struct PythonConvexHull
    {
      typedef ::opt::ConvexHull::FakeObject< T_TYPE >   t_FakeObject;
      typedef ::opt::ConvexHull::Base< t_FakeObject > t_CHBase;
      
      static void add( t_CHBase &_ch, const boost::python::tuple &_t )
      {
        __DOASSERT( len( _t ) != 3, "Convex-hull tuple must contain three objects." )
        T_TYPE object = boost::python::extract< T_TYPE >( _t[0] );
        t_FakeObject fake( object );
        ::opt::ConvexHull::NullObject::x = boost::python::extract< types::t_real >( _t[1] );
        types::t_real energy = boost::python::extract< types::t_real >( _t[2] );
        _ch.add( energy, fake );
      }
    };
 
  template<>
    struct PythonConvexHull< boost::python::object >
    {
      typedef ::opt::ConvexHull::FakeObject< boost::python::object >   t_FakeObject;
      typedef ::opt::ConvexHull::Base< t_FakeObject > t_CHBase;
      t_CHBase ch;
      
      static void add( t_CHBase &_ch, const boost::python::tuple &_t )
      {
        __DOASSERT( len( _t ) != 3, "Convex-hull tuple must contain three objects." )
        t_FakeObject fake( _t[0] );
       ::opt::ConvexHull::NullObject::x = boost::python::extract< types::t_real >( _t[1] );
        types::t_real energy = boost::python::extract< types::t_real >( _t[2] );
        _ch.add( energy, fake );
      }
    };
}

namespace opt 
{ 

  namespace ConvexHull
  {
    template<>
      std::ostream& operator<<( std::ostream &_str, const FakeObject<boost::python::object> & _o )
      {
        boost::python::str str( _o.object );
        const std::string stringobject = boost::python::extract< const std::string>(str);
        return _str << std::setw(30) << stringobject;
      }
  }
 
}

namespace PythonLaDa
{
  template< class T_TYPE > void exposeConvexHull( const std::string &_name )
  {
    using namespace boost::python;
//   const std::string fakename = "details_Fake" + _name;
//   class_< typename PythonConvexHull<T_TYPE>::t_FakeObject >
//         ( fakename.c_str(), init< typename PythonConvexHull<T_TYPE>::t_FakeObject >() );
    class_< typename PythonConvexHull<T_TYPE>::t_CHBase >( _name.c_str() )
      .def( "evaluate", &PythonConvexHull<T_TYPE>::t_CHBase::evaluate )
      .def( "add", &PythonConvexHull<T_TYPE>::add )
      .def( "__str__", &PythonConvexHull<T_TYPE>::t_CHBase::print );
  }

} //en of PythonLaDa namespace
#endif

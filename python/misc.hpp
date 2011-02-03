#ifndef __PYTHONLADA_MISC_HPP_
#define __PYTHONLADA_MISC_HPP_

#include "LaDaConfig.h"

#include <boost/python/class.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/dict.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>


#include <sstream>
#include <string>

namespace LaDa
{
  namespace Python
  {
    namespace bp = boost::python;
    template< class T_TYPE >
    std::string print( const T_TYPE &_at )
    { 
      std::ostringstream sstr;
      _at.print_out( sstr );
      return sstr.str();
    }
    template< class T_TYPE >
    std::string tostream( const T_TYPE &_at )
    { 
      std::ostringstream sstr;
      sstr << _at;
      return sstr.str();
    }

    template< class T >
      struct pickle : bp::pickle_suite
      {
        static boost::python::tuple getinitargs(T const &) { return boost::python::tuple(); }
        static bp::tuple getstate(bp::object const &_object)
        {
          T const & whatever = bp::extract<T const&>(_object);
          std::ostringstream ss;
          boost::archive::text_oarchive oa( ss );
          oa << whatever;

          return bp::make_tuple( _object.attr("__dict__"), ss.str() );
        }
        static void setstate(bp::object _out, bp::tuple state)
        {
          T & out = bp::extract<T&>(_out)();
          if (bp::len(state) == 1)
          {
            try
            {
              const std::string str = bp::extract< std::string >( state[0] );
              std::istringstream ss( str.c_str() );
              boost::archive::text_iarchive ia( ss );
              ia >> out;
            }
            catch(...)
            {
              PyErr_SetObject
              (
                PyExc_ValueError,
                ("Could not unpickle object %s\n" % (state)).ptr()
              );
              bp::throw_error_already_set();
              return;
            }
            PyErr_WarnEx( PyExc_DeprecationWarning, 
                          "This pickle is deprecated. Please repickle after unpickling.",
                          1 );
            return;
          }
          else if (bp::len(state) != 2)
          {
            PyErr_SetObject(PyExc_ValueError,
                            ("expected 2-item tuple in call to __setstate__; got %s"
                             % state).ptr()
                );
            bp::throw_error_already_set();
            return;
          }
          // restore the object's __dict__
          bp::dict d = bp::extract<bp::dict>(_out.attr("__dict__"))();
          d.update(state[0]);
          const std::string str = bp::extract< std::string >( state[1] );
          std::istringstream ss( str.c_str() );
          boost::archive::text_iarchive ia( ss );
          ia >> out;
        }
        static bool getstate_manages_dict() { return true; }
      };
  }
} // namespace LaDa
#endif // __PYTHONLADA_MISC_HPP_

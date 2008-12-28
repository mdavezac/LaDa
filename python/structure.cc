//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <sstream>

#include <boost/python.hpp>
#ifdef _MPI
# include <boost/mpi/python.hpp>
#endif
#include <boost/serialization/serialization.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <opt/types.h>
#include <opt/debug.h>
#include <mpi/macros.h>

#include "misc.hpp"
#include "xml.hpp"
#include "structure.hpp"


namespace LaDa
{
  namespace Python
  {
    namespace XML
    {
      template<> std::string nodename<Crystal::Structure>()  { return "Structure"; }
      template<> std::string nodename< Crystal::TStructure<std::string> >() 
        { return nodename<Crystal::Structure>(); }
    }

    template< class BULL >
    Crystal::Lattice& return_crystal_lattice( BULL & )
    {
      __DOASSERT( not Crystal::Structure::lattice,
                  "Lattice pointer has not been set.\n" )
      return *Crystal::Structure::lattice; 
    }
    template< class T_STRUCTURE >
      struct pickle_structure : boost::python::pickle_suite
      {
        static boost::python::tuple getinitargs( T_STRUCTURE const& _w)  
        {
          return boost::python::tuple();
        }
        static boost::python::tuple getstate(const T_STRUCTURE& _in)
        {
          std::ostringstream ss;
          boost::archive::text_oarchive oa( ss );
          oa << _in;

          return boost::python::make_tuple( ss.str() );
        }
        static void setstate( T_STRUCTURE& _out, boost::python::tuple state)
        {
          namespace bp = boost::python;
          if( bp::len( state ) != 1 )
          {
            PyErr_SetObject(PyExc_ValueError,
                            ("expected 9-item tuple in call to __setstate__; got %s"
                             % state).ptr()
                );
            bp::throw_error_already_set();
          }
          const std::string str = bp::extract< std::string >( state[0] );
          std::istringstream ss( str.c_str() );
          boost::archive::text_iarchive ia( ss );
          ia >> _out;
        }
      };

    void expose_structure()
    {
      namespace bp = boost::python;
      bp::class_< Crystal::Structure >( "Structure" )
        .def( bp::init< Crystal::Structure& >() )
        .def_readwrite( "cell",    &Crystal::Structure::cell )
        .def_readwrite( "atoms",   &Crystal::Structure::atoms )
        .def_readwrite( "k_vecs",   &Crystal::Structure::k_vecs )
        .def_readwrite( "energy",  &Crystal::Structure::energy )
        .def_readwrite( "scale",   &Crystal::Structure::scale )
        .def_readwrite( "index", &Crystal::Structure::name )
        .def( "__str__",  &print<Crystal::Structure> ) 
        .def( "fromXML",  &XML::from<Crystal::Structure> )
        .def( "toXML",  &XML::to<Crystal::Structure> )
        .def( "lattice", &return_crystal_lattice< Crystal::Structure >,
              bp::return_value_policy<bp::reference_existing_object>() )
        .def_pickle( pickle_structure< Crystal::Structure >() );
      bp::class_< Crystal::TStructure<std::string> >( "sStructure" )
        .def( bp::init< Crystal::TStructure<std::string>& >() )
        .def_readwrite( "cell",    &Crystal::TStructure<std::string>::cell )
        .def_readwrite( "atoms",   &Crystal::TStructure<std::string>::atoms )
        .def_readwrite( "energy",  &Crystal::TStructure<std::string>::energy )
        .def_readwrite( "scale",   &Crystal::TStructure<std::string>::scale )
        .def_readwrite( "index", &Crystal::TStructure<std::string>::name )
        .def( "__str__",  &print< Crystal::TStructure<std::string> > ) 
        .def( "fromXML",  &XML::from< Crystal::TStructure<std::string> > )
        .def( "toXML",  &XML::to< Crystal::TStructure<std::string> > )
        .def( "lattice", &return_crystal_lattice< Crystal::TStructure<std::string> >,
              bp::return_value_policy<bp::reference_existing_object>() )
        .def_pickle( pickle_structure< Crystal::TStructure<std::string> >() );
    }

  }
} // namespace LaDa

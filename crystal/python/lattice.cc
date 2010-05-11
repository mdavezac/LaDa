//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/object.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/errors.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <opt/types.h>
#include <opt/debug.h>
#include <crystal/structure.h>
#include <crystal/lattice.h>
#include <math/serialize.h>

#include <python/xml.hpp>
#include <python/misc.hpp>
#include "lattice.hpp"

namespace LaDa
{
  namespace Python
  {
    namespace bp = boost::python;
    namespace details
    {
      void set_as_crystal_lattice( Crystal::Lattice &_lattice )
      { 
        Crystal::Structure::lattice = &_lattice;  
        Crystal::TStructure<std::string>::lattice = &_lattice;  
      }

      template< class T_TYPE >
        void fromXML(T_TYPE &_type, const std::string &_filename )
        {
          TiXmlDocument doc( _filename ); 
          TiXmlHandle docHandle( &doc ); 
          TiXmlElement *child;
        
          __DOASSERT( not doc.LoadFile(), 
                         doc.ErrorDesc() << "\n"  
                      << "Could not load input file " << _filename  
                      << ".\nAborting.\n" ) 
          child = docHandle.FirstChild("Job").Element();
          __DOASSERT( not child,
                      "Could not find <Job> tag in " << _filename << ".\n" )
         
          child = docHandle.FirstChild("Job").FirstChild("Lattice").Element();
          __DOASSERT( not child,
                      "Could not find <Lattice> tag in " << _filename << ".\n" )

          if( child->Attribute("filename") )
          {
            const boost::filesystem::path
              n( Print::reformat_home( child->Attribute("filename") ) );
            __DOASSERT( not boost::filesystem::exists( n ),
                        n.string() + " could not be found.\n" )
            fromXML( _type, n.string() );
            return;
          }

          __DOASSERT( not _type.Load( *docHandle.FirstChild("Job").Element() ),
                         "Could not load Lattice from " + _filename + ".\n" )

          set_as_crystal_lattice( _type ); 
          _type.find_space_group();
        }

      boost::shared_ptr<Crystal::Lattice> init( std::string const &_path )
      {
        namespace bfs = boost::filesystem;
        if( not bfs::exists(_path) )
        {
          PyErr_SetString( PyExc_IOError, (_path + " does not exist.\n").c_str() );
          bp::throw_error_already_set();
          return boost::shared_ptr<Crystal::Lattice>();
        }
        if( not( bfs::is_regular(_path) or bfs::is_symlink(_path) ) )
        {
          PyErr_SetString( PyExc_IOError, (_path + " is not a file.\n").c_str() );
          bp::throw_error_already_set();
          return boost::shared_ptr<Crystal::Lattice>();
        }

        std::string text;
        opt::read_xmlfile(_path, text);
        TiXmlDocument doc;
        TiXmlHandle handle( &doc );
        doc.Parse( text.c_str() );
        TiXmlElement const *child = handle.FirstChild( "Job" ).Element();
        if( not child )
        {
          PyErr_SetString( PyExc_IOError, (_path + " does not contain Job tag.\n").c_str() );
          bp::throw_error_already_set();
          return boost::shared_ptr<Crystal::Lattice>();
        }
        child =  handle.FirstChild( "Job" ).FirstChild("Lattice").Element();
        if( not child )
        {
          PyErr_SetString( PyExc_IOError, (_path + " does not contain Lattice tag.\n").c_str() );
          bp::throw_error_already_set();
          return boost::shared_ptr<Crystal::Lattice>();
        }

        boost::shared_ptr<Crystal::Lattice> result(new Crystal::Lattice);
        if( not result->Load(*child) )
        {
          std::cerr << text << "\n";
          if( child ) std::cerr << *child <<"\n";
          PyErr_SetString( PyExc_IOError, (_path + " does not contain lattice.\n").c_str() );
          bp::throw_error_already_set();
          return boost::shared_ptr<Crystal::Lattice>();
        }

        result->find_space_group();
        set_as_crystal_lattice(*result);
        return result;
      }

    }

    math::rVector3d into_cell2(math::rVector3d const &_vec, math::rMatrix3d const &_cell)
      { return Crystal::into_cell(_vec, _cell, _cell.inverse()); }

    template< class T_STRUCTURE >
      struct pickle_lattice : bp::pickle_suite
      {
        static bp::tuple getinitargs( T_STRUCTURE const& _w)  
        {
          return bp::tuple();
        }
        static bp::tuple getstate(bp::object const &_object)
        {
          T_STRUCTURE const & structure = bp::extract<T_STRUCTURE const&>(_object);
          std::ostringstream ss;
          boost::archive::text_oarchive oa( ss );
          oa << structure;

          return bp::make_tuple( _object.attr("__dict__"), ss.str() );
        }
        static void setstate(bp::object _out, bp::tuple state)
        {
          T_STRUCTURE & out = bp::extract<T_STRUCTURE&>(_out)();
          if (bp::len(state) != 2)
          {
            PyErr_SetObject(PyExc_ValueError,
                            ("expected 2-item tuple in call to __setstate__; got %s"
                             % state).ptr()
                );
            bp::throw_error_already_set();
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

    void expose_lattice()
    {
      bp::class_< Crystal::Lattice >( "Lattice", "Defines back-bone lattice." )
        .def(bp::init< Crystal::Lattice >() )
        .def("__init__", bp::make_constructor(&details::init), "Construct lattice form xml file.\n")
        .add_property
        (
          "cell",
          make_getter(&Crystal::Lattice::cell, bp::return_value_policy<bp::return_by_value>()),
          make_setter(&Crystal::Lattice::cell, bp::return_value_policy<bp::return_by_value>()),
          "A 3x3 numpy array representing the lattice vector in cartesian units, "
          "in units of self.L{scale<lada.crystal.Lattice.scale>}."
        )
        .def_readwrite("sites", &Crystal::Lattice::sites )
        .def_readwrite("scale", &Crystal::Lattice::scale )
        .def_readwrite("space_group", &Crystal::Lattice::space_group )
        .def("__str__",  &print<Crystal::Lattice> )
        .def("fromXML",  &details::fromXML<Crystal::Lattice> )
        .def("set_as_crystal_lattice", &details::set_as_crystal_lattice )
        .def("make_primitive", &Crystal::Lattice::make_primitive,
             (bp::arg("self"), bp::arg("tolerance")=-1e0),
             "Makes lattice primitive, e.g. reduces to smallest unit-cell." )
        .def_pickle( pickle_lattice< Crystal::Lattice >() )
        .def
        ( 
          "find_space_group", 
          &Crystal::Lattice::find_space_group,
          ( bp::arg("self"), bp::arg("tolerance") = types::tolerance ),
          "Finds space-group operations (for a given tolerance), stores them in self.space_group."
        );
      bp::register_ptr_to_python< boost::shared_ptr<Crystal::Lattice> >();
      bp::def("fold_vector", &into_cell2, (bp::arg("vector"), bp::arg("cell")));
      bp::def
      (
        "fold_vector", &Crystal::into_cell,
        (bp::arg("vector"), bp::arg("cell"), bp::arg("inverse")),
        "Returns the vector folded into the given cell.\n\n"
        "@param vector: the vector to be folded.\n"
        "@type vector: numpy 3x3 float64 array.\n"
        "@param cell: the cell for which to fold.\n"
        "@type cell: numpy 3x3 float64 array.\n"
        "@param inv: the inverse of the cell for which to fold. Computed from "
        "cell if not given on input.\n"
        "@type inv: numpy 3x3 float64 array.\n"
      );
    }

  }
} // namespace LaDa

//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/errors.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

#include <opt/types.h>
#include <crystal/structure.h>
#include <crystal/lattice.h>
#include <opt/debug.h>

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

    void expose_lattice()
    {
      bp::class_< Crystal::Lattice >( "Lattice" )
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
        .def
        ( 
          "find_space_group", 
          &Crystal::Lattice::find_space_group,
          ( bp::arg("self"), bp::arg("tolerance") = types::tolerance ),
          "Finds space-group operations (for a given tolerance), stores them in self.space_group."
        );
      bp::register_ptr_to_python< boost::shared_ptr<Crystal::Lattice> >();
      bp::def( "into_cell", &Crystal::into_cell,
               (bp::arg("vector"), bp::arg("cell"), bp::arg("inverse")) );
    }

  }
} // namespace LaDa

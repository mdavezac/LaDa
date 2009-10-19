//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python.hpp>
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

    }

    void expose_lattice()
    {
      bp::class_< Crystal::Lattice >( "Lattice" )
        .def( bp::init< Crystal::Lattice >() )
        .def_readwrite( "cell",  &Crystal::Lattice::cell )
        .def_readwrite( "sites", &Crystal::Lattice::sites )
        .def_readwrite( "scale", &Crystal::Lattice::scale )
        .def_readwrite( "space_group", &Crystal::Lattice::space_group )
        .def( "__str__",  &print<Crystal::Lattice> )
        .def( "fromXML",  &details::fromXML<Crystal::Lattice> )
        .def( "set_as_crystal_lattice", &details::set_as_crystal_lattice )
        .def( "make_primitive", &Crystal::Lattice::make_primitive,
              (bp::arg("self"), bp::arg("tolerance")=-1e0),
              "Makes lattice primitive, e.g. reduces to smallest unit-cell." )
        .def
        ( 
          "find_space_group", 
          &Crystal::Lattice::find_space_group,
          ( bp::arg("self"), bp::arg("tolerance") = types::tolerance ),
          "Finds space-group operations (for a given tolerance), stores them in self.space_group."
        );
      bp::def( "into_cell", &Crystal::into_cell, (bp::arg("vector"), bp::arg("cell"), bp::arg("inverse")) );
    }

  }
} // namespace LaDa

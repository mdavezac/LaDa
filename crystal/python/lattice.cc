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
    namespace details
    {

      boost::python::list symmetry_ops( const Crystal::Lattice& _lat )
      {
        types::t_int N( _lat.space_group.point_op.getSize() );
        boost::python::list result;
        for( types::t_int i = 0; i < N; ++i )
        {
          boost::python::list inner;
          inner.append( _lat.space_group.point_op(i) );
          inner.append( _lat.space_group.trans(i) );
          result.append( inner );
        }
        return result;
      }
      
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
      using namespace boost::python;
      class_< Crystal::Lattice >( "Lattice" )
        .def( init< Crystal::Lattice >() )
        .def_readwrite( "cell",  &Crystal::Lattice::cell )
        .def_readwrite( "sites", &Crystal::Lattice::sites )
        .def_readwrite( "scale", &Crystal::Lattice::scale )
        .def( "__str__",  &print<Crystal::Lattice> )
        .def( "syms",  &details::symmetry_ops )
        .def( "fromXML",  &details::fromXML<Crystal::Lattice> )
        .def( "set_as_crystal_lattice", &details::set_as_crystal_lattice );
    // def( "StructureLattice", &return_crystal_lattice, 
    //      return_value_policy<reference_existing_object>() );
    }

  }
} // namespace LaDa

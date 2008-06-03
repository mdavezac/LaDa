//
//  Version: $Id$
//

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <iostream>
#include <algorithm>
#include <boost/lambda/lambda.hpp>

#include <opt/types.h>
#include <opt/convex_hull.h>
#include <opt/debug.h>

#include "functional_builder.h"
#include "constituent_strain.h"
#include "harmonic.h"
#include "lattice.h"
#include "structure.h"
#include "atom.h"

#include<tinyxml/tinyxml.h>

using namespace boost::python;
#include "atom.py.hpp"

BOOST_PYTHON_MODULE(atom)
{
  using namespace details;
   class_< t_Structure::t_Atoms >("VecStrings")
     .def(vector_indexing_suite< t_Site::t_Type >());

   class_< t_Atom >( "details_Atom" )
     .def( init< t_Atom >() )
     .def_readwrite( "pos",    &t_Atom::pos )
     .def_readwrite( "site",   &t_Atom::site )
     .def_readwrite( "type",   &t_Atom::type )
     .def_readwrite( "freeze", &t_Atom::freeze )
     .def( "__str__",  &print<t_Atom> ) ;
   class_< t_Site >( "details_Site" )
     .def( init< t_Site >() )
     .def_readwrite( "pos",    &t_Site::pos )
     .def_readwrite( "site",   &t_Site::site )
     .def_readwrite( "type",   &t_Site::type )
     .def_readwrite( "freeze", &t_Site::freeze )
     .def( "__str__",  &print<t_Site> ) ;


   def( "Atom", &AtomFromObject,
        return_value_policy<manage_new_object>() );
   def( "Site", &SiteFromObject,
        return_value_policy<manage_new_object>() );

   class_< t_Structure::t_Atoms >("Atoms")
     .def(vector_indexing_suite< t_Structure::t_Atoms >());
   class_< t_Structure::t_Atoms >("Sites")
     .def(vector_indexing_suite< t_Lattice::t_Sites >());

   class_< t_Structure >( "Structure" )
     .def( init< t_Structure >() )
     .def_readwrite( "cell",   &t_Structure::cell )
     .def_readwrite( "atoms",  &t_Structure::atoms )
     .def_readwrite( "energy", &t_Structure::energy )
     .def_readwrite( "scale",  &t_Structure::scale )
     .def( "__str__",  &print<t_Structure> ) 
     .def( "fromXML",  &XML::from<t_Structure> );

   class_< t_Lattice >( "Lattice" )
     .def( init< t_Lattice >() )
     .def_readwrite( "cell",  &t_Lattice::cell )
     .def_readwrite( "sites", &t_Lattice::sites )
     .def_readwrite( "scale", &t_Lattice::scale )
     .def( "__str__",  &print<t_Lattice> )
     .def( "fromXML",  &XML::from<t_Lattice> );

   typedef opt::ConvexHull::Base< PiStructure >  t_CH;
   class_< t_CH >( "ConvexHull" )
     .def( "__str__",  &t_CH::print )
     .def( "evaluate", &t_CH::evaluate )
     .def( "add",      &t_CH::add );

   class_< PiStructure >( "PiStructure" )
     .def( init< PiStructure >() )
     .def( init< types::t_int >() )
     .def( init< types::t_int, types::t_real >() )
     .def( "__str__",  &PiStructure::print )
     .def_readwrite( "index", &PiStructure::index )
     .def_readwrite( "x", &PiStructure::x );

   class_< function::Base<>::t_Container >( "details_variables" )
     .def( vector_indexing_suite< function::Base<>::t_Container >() );
     
   typedef Builder< ConstituentStrain::Harmonic::Cubic > :: t_Chemical t_Chemical;
   class_< t_Chemical >( "Chemical" )
     .def( init< t_Chemical >() )
     .def( "evaluate", &t_Chemical::evaluate );

   def( "toAtomType", &toReal );
   def( "fromAtomType", &toType );

   ExposeHarmonicRelated< ConstituentStrain::Harmonic::Cubic >();
   ExposeHarmonicRelated< ConstituentStrain::Harmonic::Tetragonal >();
}

#undef _INMODULE_

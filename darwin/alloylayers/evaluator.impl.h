//
//  Version: $Id$
//

#ifndef _DARWIN_ALLOY_LAYERS_EVALUATOR_IMPL_H_
#define _DARWIN_ALLOY_LAYERS_EVALUATOR_IMPL_H_

#include <stdexcept>       // std::runtime_error
#include <boost/lexical_cast.hpp>

#include <print/stdout.h>
#include <crystal/structure.h>
#include <crystal/atom.h>
#include <opt/fuzzy.h>
#include <opt/debug.h>
#include <mpi/mpi_object.h>

namespace GA
{
  namespace AlloyLayers
  {
#   if defined( EVALBASEHEAD ) || defined(INEVALBASE)
#     error "Macros with same names."
#   endif
#   define EVALBASEHEAD  EvaluatorBase<T_INDIVIDUAL, T_TRANSLATE, T_ASSIGN> 
#   define INEVALBASE( var ) \
       template< class T_INDIVIDUAL, \
                 template<class> class T_TRANSLATE, \
                 template<class,class> class T_ASSIGN >  var EVALBASEHEAD

    INEVALBASE( boost::shared_ptr<Crystal::Lattice> ) :: lattice( new Crystal :: Lattice );

    INEVALBASE( void ) :: init( t_Individual &_indiv )
    {
      t_Base :: init( _indiv );
      // sets structure to this object 
      translate( *current_object, structure );
    }
  
  
    INEVALBASE( bool ) :: Load( const TiXmlElement &_node )
    {
      __ASSERT( not lattice.get(), "lattice not created.\n" )
      __DOASSERT( not lattice->Load( _node ),
                  " Could not load lattice from input.\n" )
      Crystal::Structure::lattice = lattice.get();
      const TiXmlElement *parent = _node.FirstChildElement("Structure");
      __DOASSERT( not parent, "No Structure tag found in input.\n")
  
      bool loadedstructure = false;
      if ( Load_Structure( *parent ) )  return true;
  
      std::cerr << "Found attributes for constructing a layered structure...\n" 
                << "But something went wrong.\n"
                << "Continuing with standard load.\n"; 
      
      return structure.Load( *parent );
    }
  
    INEVALBASE( bool ) :: Load_Structure( const TiXmlElement &_node )
    {
      if(     ( not _node.Attribute("direction") )
          and ( not _node.Attribute("extent") )
          and ( not _node.Attribute("layersize") ) ) return false;
      if(     ( not _node.Attribute("direction") )
          or  ( not _node.Attribute("extent") ) 
          or  ( not _node.Attribute("layersize") ) ) 
      {
        std::cerr << "Either cell direction, multiplicity, "
                     "or layersize are missing on input\n";
        return false;
      }
      atat::rVector3d cdir;
      atat::rMatrix3d &cell = structure.cell;
      
      // First, Load Attributes 
      __TRYBEGIN
        std::istringstream sstr; sstr.str( _node.Attribute("direction") );
        sstr >> direction[0]; __DOASSERT( sstr.fail(), "" ) 
        sstr >> direction[1]; __DOASSERT( sstr.fail(), "" ) 
        sstr >> direction[2]; __DOASSERT( sstr.fail(), "" ) 
        __DOASSERT( atat::norm2( direction ) < types::tolerance,
                    "direction cannot be null.\n" )
        sstr.str( _node.Attribute("extent") );
        sstr >> extent[0]; __DOASSERT( sstr.fail(), "" ) 
        sstr >> extent[1]; __DOASSERT( sstr.fail(), "" ) 
        sstr >> extent[2]; __DOASSERT( sstr.fail(), "" ) 
        __DOASSERT( atat::norm2( direction ) < types::tolerance,
                    "extent cannot be null.\n" )
   
        direction = (!lattice->cell) * direction;
        
        layer_size = boost::lexical_cast<types::t_unsigned>
                                        ( _node.Attribute("layersize") );
        __DOASSERT( layer_size > 1,  "Layersize cannot be one or zero.\n" );
      __TRYEND(,"Error while loading superlattice description.\n" )
  
      // Load in scale
      if( _node.Attribute("scale") )
        _node.Attribute("scale", &structure.scale);
      else if ( structure.lattice ) 
        structure.scale = structure.lattice->scale;
      else structure.scale = 2;
  
      // Load in PI and Energy
      structure.name = "";
      if ( _node.Attribute("name") )
        structure.name = _node.Attribute("name");
      double d;
      structure.energy = 666.666;
      if( _node.Attribute("energy") )
      {
        _node.Attribute("energy", &d);
        structure.energy = (types::t_real) d;
      }
  
      bool result = Crystal :: create_epitaxial_structure( structure, direction, extent );
      __ROOTCODE
      ( 
        t_Base :: comm(),
        structure.print_xcrysden( std::cout );
        std::cout << std::endl;
      )
  
      return result;
    }
  
    INEVALBASE( inline std::string ) :: print() const
    {
      std::ostringstream sstr;
      atat::rVector3d dir = lattice->cell * direction;
      sstr << "Structure: G=" << direction 
           << ", extent=" << extent
           << ", nominal number of layers=" << layer_size << "\n"
           << structure << "\n";
      return sstr.str();
    }
  
    INEVALBASE( bool ) :: Load( t_Individual &_indiv, 
                                const TiXmlElement &_node, bool _type )
    {
      t_Object &object = _indiv.Object();
  
      if ( not object.Load( _node ) ) return false;
  
      if ( _type == GA::LOADSAVE_SHORT )
      {
        if( not _node.Attribute("string") ) return false;
        translate( _node.Attribute("string"), object );
      }
      else 
      {
        if ( not structure.Load(_node) ) return false;
        translate( structure, object );
      }
  
      assign( object, _indiv.quantities() );
  
      return true;
    }
  
    INEVALBASE( bool ) :: Save( const t_Individual &_indiv,
                                TiXmlElement &_node,
                                bool _type ) const
    {
      const t_Object &object = _indiv.Object();
      if ( _type == GA::LOADSAVE_SHORT )
      {
        std::string str;
        translate( object, str );
        _node.SetAttribute("string", str );
        return true;
      }
  
      Crystal :: Structure dummy( structure );
      translate( object, dummy );
      dummy.print_xml(_node);
  
      return true;
    }
  
#   undef EVALBASEHEAD
#   undef INEVALBASE
  
#   if defined( EVALHEAD ) || defined(INEVAL)
#     error "Macros with same names."
#   endif
#   define EVALHEAD  Evaluator<T_INDIVIDUAL, T_TRANSLATE, T_ASSIGN> 
#   define INEVAL( var ) \
       template< class T_INDIVIDUAL, \
                 template<class> class T_TRANSLATE, \
                 template<class,class> class T_ASSIGN >  var EVALHEAD
  
    INEVAL( void ) :: evaluate()
    {
      //! Calls functionals.
      Crystal::Structure copy_structure( structure );
  
      bool oldkeepdir = true;;
      if( do_dipole ) oldkeepdir = bandgap.set_keepdirectory( false );
  
      bandgap( *t_Base::current_object );
  
      if( do_dipole )
      {
        edipole( bandgap.BandGap().BandGap(), structure, *t_Base::current_object );
        bandgap.set_keepdirectory( oldkeepdir );
        if( not oldkeepdir ) bandgap.destroy_directory();
      }
  
      structure = copy_structure;
  
      // assign quantities.
      assign( *t_Base::current_object, t_Base::current_individual->quantities() );
    }
  
    INEVAL(bool) :: Load( const TiXmlElement &_node )
    {
      if( not t_Base :: Load( _node ) ) return false;
      if( not bandgap.Load( _node ) ) return false;
      edipole.Load( _node );
      return true;
    }

#   undef EVALHEAD
#   undef INEVAL
  } // namespace Layered
} // namespace GA


#endif

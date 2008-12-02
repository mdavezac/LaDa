//
//  Version: $Id$
//

#ifndef _DARWIN_GROUNDSTATES_EVALUATORBASE_IMPL_H_
#define _DARWIN_GROUNDSTATES_EVALUATORBASE_IMPL_H_

#include <stdexcept>       // std::runtime_error
#include <boost/lexical_cast.hpp>
#include <boost/bind.hpp>

#include <print/stdout.h>
#include <crystal/structure.h>
#include <crystal/atom.h>
#include <opt/fuzzy.h>
#include <opt/debug.h>
#include <mpi/mpi_object.h>
#include "../call_init.h"

namespace LaDa
{
  namespace GA
  {
    namespace GroundStates
    {
#     if defined( EVALBASEHEAD ) || defined(INEVALBASE)
#       error "Macros with same names."
#     endif
#     define EVALBASEHEAD  EvaluatorBase<T_INDIVIDUAL, T_TRANSLATE, T_ASSIGN> 
#     define INEVALBASE( var ) \
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
    
        if( not structure.Load( _node ) ) return false;

        //! now calls init member of policies.
        CallInit< t_Translate, t_This >::call( *this );
        CallInit< t_Assign, t_This >::call( *this );
 
        return true;
      }
    
      INEVALBASE( inline std::string ) :: print() const
      {
        std::ostringstream sstr;
        sstr << "Structure: \n" << structure << "\n";
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
        if ( not object.Save( _node ) ) return false;
        if ( _type == GA::LOADSAVE_SHORT )
        {
          std::string str;
          translate( object, str );
          _node.SetAttribute("string", str );
          return true;
        }
    
        translate( object, structure );
        structure.print_xml(_node);
    
        return true;
      }
    
      INEVALBASE( bool ) :: initialize( t_Individual &_indiv )
      {
        foreach( Crystal::Structure::t_Atom &atom, structure.atoms )
          atom.type = eo::rng.flip() ? 1e0: -1e0;
        translate( structure, _indiv.Object() );
        _indiv.invalidate();
        return true;
      }
    
#     undef EVALBASEHEAD
#     undef INEVALBASE
    
    } // namespace GroundStates
  } // namespace GA
} // namespace LaDa


#endif

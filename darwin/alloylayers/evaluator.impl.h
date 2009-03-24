//
//  Version: $Id$
//

#ifndef _DARWIN_ALLOY_LAYERS_EVALUATOR_IMPL_H_
#define _DARWIN_ALLOY_LAYERS_EVALUATOR_IMPL_H_

#include <stdexcept>       // std::runtime_error
#include <boost/lexical_cast.hpp>
#include <boost/bind.hpp>

#include <print/stdout.h>
#include <crystal/structure.h>
#include <crystal/atom.h>
#include <opt/tinyxml.h>
#include <opt/fuzzy.h>
#include <opt/debug.h>
#include <mpi/mpi_object.h>
#include <crystal/layerdepth.h>
#include "../call_init.h"

namespace LaDa
{
  namespace GA
  {
    namespace AlloyLayers
    {
#     if defined( EVALBASEHEAD ) || defined(INEVALBASE)
#       error "Macros with same names."
#     endif
#     define EVALBASEHEAD  EvaluatorBase<T_INDIVIDUAL, T_TRANSLATE, T_ASSIGN> 
#     define INEVALBASE( var ) \
         template< class T_INDIVIDUAL, \
                   template<class> class T_TRANSLATE, \
                   template<class,class> class T_ASSIGN >  var EVALBASEHEAD
 
      INEVALBASE( boost::shared_ptr<Crystal::Lattice> )
        :: lattice( new Crystal :: Lattice );
 
      INEVALBASE( void ) :: init( t_Individual &_indiv )
      {
        t_Base :: init( _indiv );
        // sets structure to this object 
        translate( *current_object, structure );
      }
    
    
      INEVALBASE( bool ) :: Load( const TiXmlElement &_node )
      {
        __TRYBEGIN
          boost::shared_ptr<Crystal::Lattice>
            dummy( Crystal::read_lattice( _node ) );
          lattice.swap( dummy );
        __TRYEND(, "Could not read lattice from input.\n" )
        Crystal::Structure::lattice = lattice.get();
        opt::read_tag( structure, _node, "Structure" );
        std::sort( structure.atoms.begin(), structure.atoms.end(), 
                   Crystal::LayerDepth( structure.cell ) );
 
        direction = structure.cell.get_column(0);
 
        //! now calls init member of policies.
        CallInit< t_Translate, t_This >::call( *this );
        CallInit< t_Assign, t_This >::call( *this );
 
        return true;
      }
    
      INEVALBASE( inline std::string ) :: print() const
      {
        std::ostringstream sstr;
        atat::rVector3d dir = lattice->cell * direction;
        sstr << "Structure: G=" << direction  << "\n"
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
    
#     if defined( EVALHEAD ) || defined(INEVAL)
#       error "Macros with same names."
#     endif
#     define EVALHEAD  Evaluator<T_INDIVIDUAL, T_TRANSLATE, T_ASSIGN> 
#     define INEVAL( var ) \
         template< class T_INDIVIDUAL, \
                   template<class> class T_TRANSLATE, \
                   template<class,class> class T_ASSIGN >  var EVALHEAD
    
      INEVAL( void ) :: evaluate()
      {
        //! Calls functionals.
        Crystal::Structure copy_structure( structure );
    
        bool oldkeepdir = true;;
        if( do_dipole )  oldkeepdir = bandgap.set_keepdirectory( true );
    
        bandgap( *t_Base::current_object );
 
    
        if( do_dipole )
        {
          edipole( bandgap.BandGap().BandGap(), structure, *t_Base::current_object );
          bandgap.set_keepdirectory( oldkeepdir );
          if( not oldkeepdir ) bandgap.destroy_directory();
        }
        if( do_emass )
          t_Base::current_object->emass
             =  effmass( emass, t_Base::current_object->cbm, copy_structure.cell );
        if( do_emass )
          t_Base::current_object->hmass 
             = -effmass( hmass, t_Base::current_object->vbm, copy_structure.cell );
    
        // gets concentration.
        t_Base::current_object->x = Crystal::concentration( structure, 0 );
        t_Base::current_object->y = Crystal::concentration( structure, 1 );

        structure = copy_structure;
    
        // assign quantities.
        assign( *t_Base::current_object, t_Base::current_individual->quantities() );
      }
    
      INEVAL(types::t_real) :: effmass( const Pescan::eMass& _functor,
                                        types::t_real _ref,
                                        const atat::rMatrix3d &_ocell ) const
      {
        types::t_real result = 0e0;
        Pescan::eMass::t_Output masses;
        emass( bandgap.BandGap().BandGap(), _ocell, structure, _ref, masses );
        foreach( Pescan::eMass::t_Output::value_type pair, masses )
          result += pair.second;
        return result / types::t_real( masses.size() );
      }

      INEVAL(bool) :: Load( const TiXmlElement &_node )
      {
        if( not t_Base :: Load( _node ) ) return false;
        if( not bandgap.Load( _node ) ) return false;
        if( do_dipole ) edipole.Load( _node ); // Okay if not found. Just use default values.
//       if( do_emass ) emass.Load( _node ); 
//       if( do_hmass ) hmass.Load( _node ); 
        return true;
      }

#     undef EVALHEAD
#     undef INEVAL
    } // namespace Layered
  } // namespace GA
} // namespace LaDa


#endif

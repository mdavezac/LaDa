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
  
      if( not structure.Load( _node ) ) return false;
      direction = structure.cell.get_column(0);

      __ROOTCODE
      ( 
        t_Base :: comm(),
        std::cout << structure.print_xyz( std::cout ) << std::endl;
      )

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
      std::cout << "Loading base" << std::endl;
      if( not t_Base :: Load( _node ) ) return false;
      std::cout << "Loading bandgap" << std::endl;
      if( not bandgap.Load( _node ) ) return false;
      std::cout << "Loading edipole" << std::endl;
      if( do_dipole and ( not  edipole.Load( _node ) ) ) return false;
      std::cout << "done loading" << std::endl;
      return true;
    }

    INEVAL( struct ) :: TranslateStoO
    {
      const EVALHEAD &this_;
      TranslateStoO( const EVALHEAD &_this ) : this_(_this) {}
      TranslateStoO( const TranslateStoO &_c ) : this_(_c.this_) {}
      typedef typename t_Individual :: t_IndivTraits :: t_Object t_Object;
      void operator()( const Crystal :: Structure &_s, t_Object & _o )
        { return this_.translate( _s, _o ); }
    };

    INEVAL( bool ) :: initialize( t_Individual &_indiv )
    {

      ::GA::AlloyLayers::initialize( _indiv, structure, TranslateStoO( *this ) );
      _indiv.invalidate();
      return true;
    }
#   undef EVALHEAD
#   undef INEVAL
  } // namespace Layered
} // namespace GA


#endif

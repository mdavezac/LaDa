//
//  Version: $Id$
//
#ifndef _LADA_MINIMIZER_ANY_H_
#define _LADA_MINIMIZER_ANY_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/variant.hpp>
#include <boost/variant/apply_visitor.hpp>
#include <boost/mpl/vector.hpp>

#include <opt/tinyxml.h>

#include "frprmn.h"
#include "gsl_mins.h"
#include "decoupled.h"

namespace LaDa
{
  namespace Minimizer
  {
    class Decoupled;

    //! Wraps around other minimimizers to provided single Load/Launch.
    class Any
    {
      public:
        //! Constructor.
        Any() {}
        //! Copy Constructor.
        Any( const Any& _c ) : minimizer( _c.minimizer ) {}


        //! Calls minimizer.
        template< class T_FUNCTOR >
          typename T_FUNCTOR :: t_Return 
            operator()( const T_FUNCTOR& _a, typename T_FUNCTOR :: t_Arg& _b );
        
        //! Load Frpr or Gsl minimizer from XML.
        bool Load( const TiXmlElement& _node );

      protected:
        //! Vector Types.
        typedef boost::mpl::vector< Frpr, Gsl, Decoupled > t_Minimizers;
        //! A points to a Gsl minimizer.
        boost::make_variant_over< t_Minimizers > minimizer;
    };

    namespace details
    {
      template< class T_FUNCTION >
        struct OpVisitor : public boost::static_visitor< typename T_FUNCTION :: t_Return > 
        {
          T_FUNCTION &func_;
          typename T_FUNCTION :: t_Arg &arg_;
        
          OpVisitor   ( const T_FUNCTION &_func, typename T_FUNCTION :: t_Arg &_arg )
                    : func_( _func ), arg_( _arg ) {}
          OpVisitor( const OpVisitor &_c ) : func_( _c.func_ ), arg_( _c.arg_ ) {}

          template< class T_MINIMIZER >
            typename T_FUNCTION :: t_Return operator()( const T_MINIMIZER& _minimizer ) const
              { return _minimizer( func_, arg_ ); }
        };
    }

    template< class T_FUNCTOR >
      typename T_FUNCTOR :: t_Return 
        Any :: operator()( const T_FUNCTOR& _a, typename T_FUNCTOR :: t_Arg& _b )
        {
          typedef typename T_FUNCTOR :: t_Return t_Return;
          typedef typename T_FUNCTOR :: t_Arg t_Arg;
          return boost::apply_visitor( details::OpVisitor<T_FUNCTOR>( _a, _b ), minimizer );
        }
  }
} 
#endif

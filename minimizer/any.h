//
//  Version: $Id$
//
#ifndef _LADA_MINIMIZER_ANY_H_
#define _LADA_MINIMIZER_ANY_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/function.hpp>
#include <boost/type_traits/remove_pointer.hpp> 

#include <opt/tinyxml.h>

#include "frprmn.h"
#include "gsl_mins.h"

namespace LaDa
{
  namespace Minimizer
  {
    //! Wraps around other minimimizers to provided single Load/Launch.
    class Any
    {
      public:
        //! Constructor.
        Any() {}
        //! Copy Constructor.
        Any( const Any& _c ) : gsl_( _c.gsl_ ), frpr_( _c.frpr_ ) {}


        template< class T_FUNCTOR >
          typename T_FUNCTOR :: t_Return 
            operator()( const T_FUNCTOR& _a, typename T_FUNCTOR :: t_Arg& _b )
            { 
              if( bool( gsl_ ) ) return (*gsl_)( _a, _b ); 
              if( bool( frpr_ ) ) return (*frpr_)( _a, _b ); 
              __ASSERT( true, "No Minimizer is present.\n" )
            }
        
        //! Load Frpr or Gsl minimizer from XML.
        bool Load( const TiXmlElement& _node );

      protected:
        //! A points to a Gsl minimizer.
        boost::shared_ptr< Gsl > gsl_;
        //! A points to a Frpr minimizer.
        boost::shared_ptr< Frpr > frpr_;
    };

  }
} 
#endif

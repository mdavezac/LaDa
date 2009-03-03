//
//  Version: $Id$
//

#ifndef _LADA_PYTHON_FUNCTION_HPP_
#define _LADA_PYTHON_FUNCTION_HPP_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/python/class.hpp>
#include <boost/python/object.hpp>
#include <boost/python/list.hpp>
#include <boost/python/extract.hpp>
#include <boost/type_traits/remove_pointer.hpp>

#include <opt/debug.h>

namespace LaDa
{
  namespace Python
  {
    class Function
    {
      public:
        //! Type of the return.
        typedef types :: t_real t_Return;
        //! Type of the container.
        typedef std::vector< types :: t_real > t_Arg;
        //! Type of the gradient argument.
        typedef types :: t_real* t_GradientArg;
        //! Constructor.
        Function( const boost :: python :: object &_object ) : object_( _object ) {}
        t_Return operator()( const t_Arg& _arg ) const
        {
          const boost::python::object result( object_( _arg ) );
          return boost::python::extract< t_Return>( result  ); 
        }
        //! calls gradient.
        t_Return gradient( const t_Arg& _arg, t_GradientArg _gradient ) const
        {
          namespace bp = boost :: python;

          std::vector< t_Arg::value_type > gradient( _arg.size(), 0 );
          // calls gradient function.
          object_.attr("gradient")( _arg, gradient );
          std::vector< t_Arg::value_type > :: iterator i_var = gradient.begin();
          std::vector< t_Arg::value_type > :: iterator i_var_end = gradient.end();
          for(; i_var != i_var_end; ++i_var, ++_gradient ) *_gradient += *i_var;
        }

      protected:
        //! References an object.
        const boost :: python :: object& object_;
    };
  }
} // namespace LaDa
#endif

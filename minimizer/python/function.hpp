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
        //! callable operator.
        t_Return operator()( const t_Arg& _arg ) const
        {
          const boost::python::object result( object_( _arg ) );
          return boost::python::extract< t_Return>( result  ); }
        //! calls gradient.
        t_Return gradient( const t_Arg& _arg, t_GradientArg _gradient ) const
        {
          namespace bp = boost :: python;
          // copy to gradient
          bp::list grad;
          t_GradientArg i_grad( _gradient );
          const size_t n(_arg.size() ); 
          for( size_t i(0); i < n; ++i, ++i_grad ) grad.append( *i_grad );
          // calls gradient function.
          object_.attr("gradient")( _arg, grad );
          // copies gradient back
          i_grad = _gradient;
          for( size_t i(0); i < n; ++i, ++i_grad )
            *i_grad = bp::extract< boost::remove_pointer< t_GradientArg > :: type >( grad[i] );
        }

      protected:
        //! References an object.
        const boost :: python :: object& object_;
    };
  }
} // namespace LaDa
#endif

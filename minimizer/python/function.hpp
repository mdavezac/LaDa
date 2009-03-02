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
        //! callable operator.
        t_Return operator()( const boost::python::list& _arg ) const
        {
          const boost::python::object result( object_( _arg ) );
          return boost::python::extract< t_Return>( result  ); 
        }
        //! callable operator.
        t_Return operator()( const t_Arg& _arg ) const
        {
          boost::python::list arg;
          t_Arg :: const_iterator i_var = _arg.begin();
          t_Arg :: const_iterator i_var_end = _arg.end();
          for(; i_var != i_var_end; ++i_var ) arg.append( *i_var );
          return (*this)( arg );
        }
        //! calls gradient.
        t_Return gradient( const boost::python::list& _arg, boost::python::list &_gradient ) const
        {
          namespace bp = boost :: python;

          // calls gradient function.
          object_.attr("gradient")( _arg, _gradient );
        }
        //! calls gradient.
        t_Return gradient( const t_Arg& _arg, t_GradientArg _gradient ) const
        {
          namespace bp = boost :: python;
          boost::python::list arg, grad;
          t_Arg :: const_iterator i_var = _arg.begin();
          t_Arg :: const_iterator i_var_end = _arg.end();
          t_GradientArg i_grad( _gradient );
          for(; i_var != i_var_end; ++i_var, ++i_grad )
          {
            arg.append( *i_var );
            grad.append( *i_grad );
          }

          gradient( arg, grad );
          // copies gradient back
          const size_t n( _arg.size() ); 
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

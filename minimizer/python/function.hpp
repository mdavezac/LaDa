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
#include <boost/python/errors.hpp>

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
          PyObject* error = PyErr_Occurred();
          if( error )
          {
            std::cerr << "Python raised an exception. Could not obtain gradient.\n";
            boost::python::throw_error_already_set();
            return 0;
          }
          __TRYBEGIN
            const t_Return value = boost::python::extract< t_Return>( result  ); 
            return value;
          __TRYEND(,"Could not extract real value from return of python object.\n")
        }
        //! calls gradient.
        void gradient( const t_Arg& _arg, t_GradientArg _gradient ) const
        {
          __TRYBEGIN
          namespace bp = boost :: python;

          typedef std::vector< t_Arg::value_type > t_vector;
          t_vector gradient( _arg.size(), 0);
          // calls gradient function.
          const boost::python::object result( object_.attr("gradient")( _arg, gradient ) ); 
          PyObject* error = PyErr_Occurred();
          if( error )
          {
            std::cerr << "Python raised an exception. Could not obtain gradient.\n";
            bp::throw_error_already_set();
          }
          __TRYBEGIN
            gradient = bp::extract< t_vector >( result );
          __TRYEND(,"Could not extract lada.opt.cReal from return of python object.\n")
          std::vector< t_Arg::value_type > :: iterator i_var = gradient.begin();
          std::vector< t_Arg::value_type > :: iterator i_var_end = gradient.end();
          for(; i_var != i_var_end; ++i_var, ++_gradient ) *_gradient += *i_var;
          __TRYEND(, "Error in gradient.\n" )
        }

      protected:
        //! References an object.
        const boost :: python :: object& object_;
    };
  }
} // namespace LaDa
#endif

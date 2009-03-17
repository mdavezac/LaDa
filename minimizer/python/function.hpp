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
          __TRYBEGIN
          namespace bp = boost :: python;

          typedef std::vector< t_Arg::value_type > t_vector;
          t_vector gradient( _arg.size(), 0);
          // calls gradient function.
          gradient = bp::extract< t_vector >( object_.attr("gradient")( _arg, gradient ) );
          std::vector< t_Arg::value_type > :: iterator i_var = gradient.begin();
          std::vector< t_Arg::value_type > :: iterator i_var_end = gradient.end();
          std::cout << "python ";
          t_GradientArg grad( _gradient );
          for(; i_var != i_var_end; ++i_var, ++grad )
          { 
            *grad += *i_var;
            std::cout << *i_var << " ";
          }
          std::cout << "\n";
          __TRYEND(, "Error in gradient.\n" )
        }

      protected:
        //! References an object.
        const boost :: python :: object& object_;
    };
  }
} // namespace LaDa
#endif

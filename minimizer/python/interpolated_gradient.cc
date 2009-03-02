//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python/def.hpp>
#include <boost/python/object.hpp>
#include <boost/python/list.hpp>
#include <boost/python/extract.hpp>

#include "../cgs.h"
#include "../interpolated_gradient.h"
#include "function.hpp"


namespace LaDa
{
  namespace Python
  {
    void interpolated_gradient( const boost::python::object &_function, 
                                const boost::python::list &_arg, 
                                boost::python::list &_gradient, 
                                const size_t _n, 
                                const Function :: t_Return _stepsize,
                                const size_t _itermax,
                                const types::t_real _tolerance,
                                const bool _verbose )
    {
      namespace bp = boost::python;
      typedef Fitting :: Cgs t_Cgs;
      typedef Function t_Function;
      const t_Function function( _function );
      t_Cgs cgs;
      cgs.verbose = _verbose;
      cgs.itermax = _itermax;
      cgs.tolerance = _tolerance;
      typedef boost::remove_pointer<t_Function :: t_GradientArg > :: type t_Type;
      std::vector< t_Type > gradient( boost::python::len( _arg ), 0 );
      t_Function :: t_Arg arg( boost::python::len( _arg ), 0 );
      t_Function :: t_Arg :: iterator i_var = arg.begin();
      t_Function :: t_Arg :: iterator i_var_end = arg.end();
      for(size_t i(0); i_var != i_var_end; ++i_var, ++i )
        *i_var = bp :: extract< t_Type >( _arg[i] );

      __DOASSERT( bp::len( _gradient ) != bp::len( _arg ),
                  "Gradient and argument have different sizes.\n" )
      
      Minimizer :: interpolated_gradient< t_Function, t_Cgs >
      (
        function, 
        arg,
        cgs,
        &(gradient[0]),
        _n,
        _stepsize
      );
      i_var = gradient.begin();
      i_var_end = gradient.end();
      for(size_t i(0); i_var != i_var_end; ++i_var, ++i ) 
        _gradient[i] = bp::extract< t_Type >( _gradient[i] ) + (*i_var);
    }

    void expose_interpolated_gradient()
    {
      namespace bp = boost::python;
      bp::def
      ( 
        "interpolated_gradient", 
        &interpolated_gradient,
        (
          bp::arg("self"), 
          bp::arg("arg"), 
          bp::arg("gradient"), 
          bp::arg("n") = 1,
          bp::arg("stepsize") = 1e-3,
          bp::arg("itermax") = 50,
          bp::arg("tolerance") = 1e-12,
          bp::arg("verbosity") = false
        ),
        "Interpolates gradient from calls to self.\n"
        "self should be callable.\n"
        "arg should be convertible to a (c++) vector of reals.\n"
        "gradient should be a list of reals.\n"
        "n is the order of the interpolation.\n"
        "stepsize are the steps taken for the interpolation.\n"
        "itermax is the maximum number of iterations when performing the fit.\n"
        "tolerance is the tolerance of the fit.\n"
        "If verbosity is true, outputs fitting stuff.\n"
      );
    }
  } // namespace Python
} // namespace LaDa

//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python/def.hpp>

#include "../cgs.h"
#include "../interpolated_gradient.h"
#include "function.hpp"


namespace LaDa
{
  namespace Python
  {
    void interpolated_gradient( const Function &_function, 
                                const Function :: t_Arg &_arg, 
                                boost :: python :: list &_gradient, 
                                const size_t _n, 
                                const Function :: t_Return _stepsize,
                                const size_t _itermax,
                                const types::t_real _tolerance,
                                const bool _verbose )
    {
      typedef Fitting :: Cgs t_Cgs;
      typedef Function t_Function;
      t_Cgs cgs;
      cgs.verbose = _verbose;
      cgs.itermax = _itermax;
      cgs.tolerance = _tolerance;
      typedef boost::remove_pointer<t_Function :: t_GradientArg > :: type t_Type;
      std::vector< t_Type > gradient( _arg.size(), 0 );
      
      Minimizer :: interpolated_gradient< t_Function, t_Cgs >
      (
        _function, 
        _arg,
        cgs,
        &(gradient[0]),
        _n,
        _stepsize
      );
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

#include "LaDaConfig.h"

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
                                const std::vector<types::t_real> &_arg, 
                                std::vector<types::t_real> &_gradient, 
                                const size_t _n, 
                                const Function :: t_Return _stepsize,
                                const size_t _itermax,
                                const types::t_real _tolerance,
                                const bool _verbose )
    {
      __DOASSERT( _gradient.size() != _arg.size(), "Argurments and gradients have different sizes.\n" )
      namespace bp = boost::python;
      typedef Fitting :: Cgs t_Cgs;
      typedef Function t_Function;
      const t_Function function( _function );
      t_Cgs cgs;
      cgs.verbose = _verbose;
      cgs.itermax = _itermax;
      cgs.tolerance = _tolerance;
      
      Minimizer :: interpolated_gradient< t_Function, t_Cgs >
      (
        function, 
        _arg,
        cgs,
        &(_gradient[0]),
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
          bp::arg("verbose") = false
        ),
        "Interpolates gradient using numerical derivatives.\n"
        "@param self: callable object taking arg on input and returning a scalar.\n"
        "@param arg: a vector of real numbers.\n"
        "@type arg: c++ vector\n"
        "@param gradient: list of gradients.\n"
        "@type gradient: c++ vector\n"
        "@param n:  order of the interpolation.\n"
        "@type n: integer\n"
        "@param stepsize: size of the steps taken for the interpolation (default: 1e-3).\n"
        "@type stepsize: float\n"
        "@param itermax:  maximum number of iterations when performing the fit (default: 50).\n"
        "@type itermax: integer\n"
        "@param tolerance: tolerance of the fit (default: 1e-12).\n"
        "@type tolerance: float\n"
        "@param verbose: True gets you more output (default: False).\n"
        "@type verbose: Boolean\n"
      );
    }
  } // namespace Python
} // namespace LaDa

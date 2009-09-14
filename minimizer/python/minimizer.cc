//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif


#include <boost/mpl/vector.hpp>


#include "minimizer.hpp"


namespace LaDa
{
  namespace Python
  {

    void Minimizer :: load_()
    {
      TiXmlElement fakexml( "Minimizer" );
      fakexml.SetAttribute( "type", type_ );
      fakexml.SetDoubleAttribute( "tolerance", tolerance_ );
      fakexml.SetAttribute( "itermax", itermax_ );
      fakexml.SetDoubleAttribute( "linetolerance", linetolerance_ );
      fakexml.SetDoubleAttribute( "linestep", linestep_ );
      fakexml.SetAttribute( "strategy", strategy_ );
      fakexml.SetAttribute( "verbose", verbose_ );
      fakexml.SetAttribute( "uncertainties", uncertainties_ );
      fakexml.SetAttribute( "up", up_ );
      fakexml.SetAttribute( "gradient", use_gradient_ ? "true": "false" );
      __DOASSERT( not minimizer_.Load( fakexml ), "Could not load minimizer " << type_ << ".\n" )
    }

    void expose_minimizer()
    {
      namespace bp = boost::python;
      bp::class_< Minimizer >
      ( 
         "Minimizer",
         "Will minimizer a function object. The function object must be callable,\n"
         "and must define a gradient function as well.\n"
         "The minimizer can be chosen from any of gsl, original vff, and minuit2 minimizers.\n"
      )
#       define __GETSET__(NAME, DOC) \
          .add_property( #NAME, &Minimizer::get_ ## NAME, &Minimizer::set_ ## NAME, DOC )
        __GETSET__( type, "Name of the minimizer to use." )
        __GETSET__( tolerance, "Convergence criteria." )
        __GETSET__( itermax, "Maximum number of iterations." )
        __GETSET__( linetolerance, "Line tolerance for conjugate-gradient methods." )
        __GETSET__( linestep, "Line steo for conjugate-gradient methods." )
        __GETSET__( strategy, "slowest/slow/fast strategies for Minuit2 minimizers." )
        __GETSET__( verbose, "If true, will print verbose output." )
        __GETSET__( use_gradient, "If true will use gradient from object, if false, minuit2 will compute gradient." )
#       undef __GETSET__
        .def( "__call__", &Minimizer::operator() )
        .def(
              "set", &Minimizer::set, 
              (
                bp::arg("type") = "gsl_bfgs2",
                bp::arg("convergence") = 1e-6,
                bp::arg("itermax") = 50,
                bp::arg("linetolerance") = 1e-2,
                bp::arg("linestep") = 1e-1,
                bp::arg("strategy") = "slow",
                bp::arg("verbose") = false,
                bp::arg("uncertainties") = 0.1,
                bp::arg("up") = 1,
                bp::arg("gradient") = true
              ),
              "Sets parameters for the optimizers. "
              "Not all parameters are needed by all optimizers."
            );
    }
  } // namespace Python
} // namespace LaDa
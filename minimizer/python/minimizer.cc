#include "LaDaConfig.h"

#include<boost/python/class.hpp>
#include<boost/python/tuple.hpp>
#include<boost/python/make_constructor.hpp>

#include <boost/mpl/vector.hpp>

#include <python/misc.hpp>
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
      LADA_DO_NASSERT( not minimizer_.Load( fakexml ), "Could not load minimizer " << type_ << ".\n" )
    }

    t_Minimizer* make_internal(Minimizer const &_in)
    { 
      t_Minimizer *result(new t_Minimizer(_in.minimizer_));
      if(result == NULL)
      {
        PyErr_SetString(PyExc_RuntimeError, "Could not create internal C++ minimizer.");
        bp::throw_error_already_set();
      }
      return result;
    }

    boost::python::tuple expose_minimizer()
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
        __GETSET__( linestep, "Line step for conjugate-gradient methods." )
        __GETSET__( strategy, "slowest/slow/fast strategies for Minuit2 minimizers." )
        __GETSET__( up, "Not sure. For Minuit2 minimizers only." )
        __GETSET__( uncertainties, "Not sure. For Minuit2 minimizers only, and "
                                   "zeps of frpmn minimizer." )
        __GETSET__( verbose, "If true, will print verbose output." )
        __GETSET__( use_gradient, "If true will use gradient from object, if false, minuit2 will compute gradient." )
#       undef __GETSET__
        .def( "__call__", &Minimizer::operator() )
        .def(
              "set", &Minimizer::set, 
              (
                bp::arg("type") = "gsl_bfgs2",
                bp::arg("tolerance") = 1e-6,
                bp::arg("itermax") = 50,
                bp::arg("linetolerance") = 1e-2,
                bp::arg("linestep") = 1e-1,
                bp::arg("strategy") = "slow",
                bp::arg("verbose") = false,
                bp::arg("uncertainties") = 0.1,
                bp::arg("up") = 1,
                bp::arg("use_gradient") = true
              ),
              "Sets parameters for the optimizers. "
              "Not all parameters are needed by all optimizers."
            )
        .def("_copy_to_cpp", &Minimizer::copy_to_internal)
        .def_pickle( Python::pickle<Minimizer>() );

      bp::class_<t_Minimizer>("_CppMinimizer", bp::no_init)
        .def("__init__", bp::make_constructor(&make_internal));

      return bp::make_tuple
        (
          "frprmn"
#       ifdef LADA_WITH_GSL
          , "gsl_bfgs", "gsl_bfgs2", "gsl_fr", "gsl_pr", "gsl_sd"
#       endif
#       ifdef LADA_WITH_MINUIT2
          , "minuit2"
#       endif
        );
    }
  } // namespace Python
} // namespace LaDa

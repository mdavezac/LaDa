//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif


#include <boost/mpl/vector.hpp>

#include <opt/debug.h>

#include "../minuit2.h"
#include "../frprmn.h"
#include "../gsl_mins.h"
#include "function.hpp"


namespace LaDa
{
  namespace Python
  {
    class Minimizer
    {
      public:
        Minimizer() : type_( "gsl_bfgs2" ),
                      tolerance_( 1e-6 ),
                      itermax_( 50 ),
                      linetolerance_( 1e-2 ),
                      linestep_( 1e-3 ),
                      strategy_( "fast" ),
                      verbose_( false ) {}

#       define __GETSET__( TYPE, NAME ) \
          TYPE get_ ## NAME () const { return NAME ## _; }\
          void set_ ## NAME( const TYPE& _ ## NAME ) \
          {\
            NAME ## _ = _ ## NAME; \
            __TRYBEGIN \
            load_(); \
            __TRYEND(,"Could not set parameter " #NAME ".\n") \
          }

        __GETSET__( std::string, type)
        __GETSET__( types::t_real, tolerance )
        __GETSET__( size_t, itermax )
        __GETSET__( types::t_real, linetolerance )
        __GETSET__( types::t_real, linestep )
        __GETSET__( std::string, strategy )
        __GETSET__( bool, verbose )
#       undef __GETSET__

        Function :: t_Return operator()( const boost::python::object &_function,
                                         Function :: t_Arg& _arg ) const 
        {
          const Function function( _function );
          return minimizer_( function, _arg ); 
        }
        void set( const std::string& _type, 
                  types::t_real _tolerance, 
                  size_t _itermax, 
                  types::t_real _linetolerance,
                  types::t_real _linestep,
                  const std::string& _strategy,
                  bool _verbose )
        {
          type_ = _type;
          tolerance_ = _tolerance;
          itermax_ = _itermax;
          linetolerance_ = _linetolerance;
          linestep_ = _linestep;
          strategy_ = _strategy;
          verbose_ = _verbose;
          load_();
        }

      protected:
        void load_();
        std::string type_;
        types::t_real tolerance_;
        size_t itermax_;
        types::t_real linetolerance_;
        types::t_real linestep_;
        std::string strategy_;
        bool verbose_;
        typedef LaDa::Minimizer::Variant
                < 
                  boost::mpl::vector
                  <
                    LaDa::Minimizer::Frpr, 
                    LaDa::Minimizer::Gsl, 
                    LaDa::Minimizer::Minuit2
                  > 
                > t_Minimizer;
        t_Minimizer minimizer_;
    };

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
                bp::arg("verbose") = false
              ),
              "Sets parameters for the optimizers. "
              "Not all parameters are needed by all optimizers."
            );
    }
  } // namespace Python
} // namespace LaDa

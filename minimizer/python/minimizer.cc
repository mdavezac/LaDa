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

#       define __getset__ ( TYPE, NAME ) \
          TYPE get_ ## NAME () const { return NAME ## _; }\
          void set_ ## NAME( const TYPE& _ ## NAME ) \
          {\
            _ ## NAME = NAME ## _;\
            __TRYBEGIN \
            load_(); \
            __TRYEND(,"Could not set parameter " ## NAME ## " ") \
          }

        __getset__( std::string, type)
        __getset__( types::t_real, tolerance )
        __getset__( size_t, itermax )
        __getset__( types::t_real, linetolerance )
        __getset__( types::t_real, linestep )
        __getset__( std::string, strategy )
#       undef __getset__

        Function :: t_Return operator()( const Function& _function,
                                         const Function :: t_Arg& _arg ) const 
          { return minimizer_( _function, _arg ); }

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
      fakexml.SetAttribute( "type", minimizer_type );
      fakexml.SetDoubleAttribute( "tolerance", tolerance );
      fakexml.SetAttribute( "itermax", itermax );
      fakexml.SetDoubleAttribute( "linetolerance", linetol );
      fakexml.SetDoubleAttribute( "linestep", linestep );
      fakexml.SetAttribute( "strategy", strategy );
      fakexml.SetAttribute( "verbose", verbosity >= outermin ? "true": "false" );
      __DOASSERT( not minimizer_.Load( fakexml ), "Could not load minimizer " << type_ << ".\n" )
    }

    void expose_minimizer()
    {
      namespace bp = boost::python;
      bp::class_< Minimizer >( "Minimizer" )
#       define __getset__(NAME, DOC) \
          .add_property( #NAME, &Minimizer::get_ ## NAME, &Minimizer::set_ ## NAME, DOC )
        __getset__( type, "" )
        __getset__( tolerance, "" )
        __getset__( itermax, "" )
        __getset__( linetolerance, "" )
        __getset__( linestep, "" )
        __getset__( strategy, "" )
        __getset__( verbose, "" )
#       undef __getset__
        .def( "__call__", &Minimizer::operator() );
    }
  } // namespace Python
} // namespace LaDa

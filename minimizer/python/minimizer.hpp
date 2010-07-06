#ifndef _LADA_PYTHON_MINIMIZER_HPP_
#define _LADA_PYTHON_MINIMIZER_HPP_

#include "LaDaConfig.h"

#include <boost/serialization/serialization.hpp>

#include <opt/debug.h>

#include "../frprmn.h"
#if defined(LADA_WITH_MINUIT2) or defined(LADA_WITH_GSL)
# include "../variant.h"
#endif
#ifdef LADA_WITH_MINUIT2
# include "../minuit2.h"
#endif
#ifdef LADA_WITH_GSL
# include "../gsl_mins.h"
#endif
#include "function.hpp"

namespace LaDa
{
  namespace Python
  {
    void expose_minimizer();

    //! Declaration for the python minimizer.
    class Minimizer
    {
      friend class boost::serialization::access;
      public:
        Minimizer() : type_( "gsl_bfgs2" ),
                      tolerance_( 1e-6 ),
                      itermax_( 50 ),
                      linetolerance_( 1e-2 ),
                      linestep_( 1e-3 ),
                      strategy_( "fast" ),
                      verbose_( false ),
                      uncertainties_(0.1),
                      up_(1),
                      use_gradient_( true ) {}

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
        __GETSET__( types::t_real, uncertainties )
        __GETSET__( types::t_real, up )
        __GETSET__( bool, use_gradient )
#       undef __GETSET__

        Function :: t_Return operator()( const boost::python::object &_function,
                                         Function :: t_Arg& _arg ) const 
        {
          try
          {
            const Function function( _function );
            return minimizer_( function, _arg ); 
          }
          catch (...)
          {
            PyErr_SetString( PyExc_RuntimeError, "Error encountered while minimizing from C.\n");
            return -1;
          }
        }
        //! Call minimizer, without python callback.
        template< class T_FUNCTOR >
          typename T_FUNCTOR :: t_Return 
            call( const T_FUNCTOR& _a, typename T_FUNCTOR :: t_Arg& _b ) const
            {
              try
              {
                return minimizer_( _a, _b ); 
              }
              catch (...)
              {
                PyErr_SetString( PyExc_RuntimeError, "Error encountered while minimizing from C.\n");
                return -1;
              }
            }

        void set( const std::string& _type, 
                  types::t_real _tolerance, 
                  size_t _itermax, 
                  types::t_real _linetolerance,
                  types::t_real _linestep,
                  const std::string& _strategy,
                  bool _verbose,
                  types::t_real _uncertainties,
                  types::t_real _up,
                  bool _ug )
        {
          try
          {
            type_ = _type;
            tolerance_ = _tolerance;
            itermax_ = _itermax;
            linetolerance_ = _linetolerance;
            linestep_ = _linestep;
            strategy_ = _strategy;
            verbose_ = _verbose;
            uncertainties_ = _uncertainties;
            up_ = _up;
            use_gradient_ = _ug;
            load_();
          }
          catch (...)
          {
            PyErr_SetString( PyExc_RuntimeError, "Error setting minimizer parameters.\n");
          }
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
        types::t_real uncertainties_;
        types::t_real up_;
        bool use_gradient_;
#       if defined(LADA_WITH_MINUIT2) or defined(LADA_WITH_GSL)
          typedef LaDa::Minimizer::Variant
                  < 
                    boost::mpl::vector
                    <
                      LaDa::Minimizer::Frpr  
#                     ifdef LADA_WITH_GSL
                        , LaDa::Minimizer::Gsl  
#                     endif
#                     ifdef LADA_WITH_MINUIT2
                        , LaDa::Minimizer::Minuit2
#                     endif
                    > 
                  > t_Minimizer;
#       else
          typedef LaDa::Minimizer::Frpr t_Minimizer;
#       endif
        t_Minimizer minimizer_;

        //! Serializes a lattice.
        template<class ARCHIVE> void serialize(ARCHIVE & _ar, const unsigned int _version)
        {
           _ar & type_; _ar & tolerance_; _ar & itermax_; _ar & linetolerance_;
           _ar & linestep_; _ar & strategy_; _ar & verbose_; _ar & uncertainties_;
           _ar & up_; _ar & use_gradient_;
        }
    };
  }
} // namespace LaDa
#endif

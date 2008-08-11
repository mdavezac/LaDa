//
//  Version: $Id$
//
#ifndef _CE_SEPMAPPINGS_H_
#define _CE_SEPMAPPINGS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/lambda/lambda.hpp>

#include <opt/types.h>
#include <opt/debug.h>

namespace CE
{
  //! Contains mapping policies for fixed-lattice separable functions.
  namespace Mapping
  {
    //! \brief Allows different types of mapping from confs to coef parameters.
    //! \detail This mapping is equivalent to VectorPlus, eg (1,0..), (0,1,....),
    //!          and so on.
    template< size_t DIM > class VectorPlus 
    {
      public:
        //! A D dimensional mapping.
        const static size_t D = DIM;
        //! Applies function itself.
        template< class T_CONF, class T_ITCOEF, class T_OUT >
          const static void apply( const T_CONF &_conf,
                                   const T_ITCOEF &_coef, 
                                   T_OUT &_out )
            { _out *= *( _coef + typename T_ITCOEF::difference_type( _conf ) ); }
        //! Returns the normalized, coef-less vector.
        template< class T_CONF, class T_OUT >
          const static void add_tovec( const T_CONF &_conf, T_OUT &_out, 
                                       typename T_OUT::value_type _s  
                                           = typename T_OUT::value_type(1) )
          { _out[ size_t( _conf ) ] += typename T_OUT::value_type(_s); }
        //! Normalizes vector.
        template< class T_ITCOEF, class T_NORM >
          const static void apply( T_ITCOEF &_coef,  T_NORM &_out )
          {
            namespace bl = boost::lambda;
            types::t_real norm(0);
            std::for_each( _coef, _coef + D, bl::var(norm) += bl::_1 * bl::_1 );
            if( Fuzzy::is_zero( norm) ) return;
            norm = std::sqrt( norm );
            _out *= norm;
            norm = T_NORM(1) / norm;
            std::for_each( _coef, _coef + D, bl::_1 *= bl::constant(norm) );
          }
    };
    //! \brief Allows different types of mapping from confs to coef parameters.
    //! \detail This mapping is equivalent to VectorPlus, with one constant
    //!         vector, and all other vectors with a single non-zero component.
    template< size_t DIM > class VectorDiff
    {
      public:
        //! A D dimensional mapping.
        const static size_t D = DIM;
        //! Applies functions with appropriate coef.
        template< class T_CONF, class T_ITCOEF, class T_OUT >
          const static void apply( const T_CONF &_conf,
                                   const T_ITCOEF &_coef, 
                                   T_OUT &_out )
          {
            _out *= *_coef;
            if( Fuzzy::is_zero( _conf ) ) return;
            _out *= *( _coef + typename T_ITCOEF::difference_type( _conf ) );
          }
        //! Returns the normalized, coef-less vector.
        template< class T_CONF, class T_OUT >
          const static void add_tovec( const T_CONF &_conf, T_OUT &_out, 
                                       const typename T_OUT::value_type _s 
                                               = typename T_OUT::value_type(1) )
          {
            _out[0] += typename T_OUT::value_type(_s);
            if( Fuzzy::is_zero( _conf ) ) return;
            _out[ size_t( _conf ) ] +=  typename T_OUT::value_type(_s); 
          }
        //! Normalizes vector.
        template< class T_ITCOEF, class T_NORM >
          const static void apply( T_ITCOEF &_coef,  T_NORM &_out )
          {
            namespace bl = boost::lambda;
            types::t_real norm(0);
            std::for_each( _coef, _coef + D, bl::var(norm) += bl::_1 * bl::_1 );
            if( Fuzzy::is_zero( norm) ) return;
            norm = std::sqrt( norm );
            _out *= norm;
            norm = T_NORM(1) / norm;
            std::for_each( _coef, _coef + D, bl::_1 *= bl::constant(norm) );
          }
    };
  } // end of Mapping namespace.
} // end of CE namespace.
#endif

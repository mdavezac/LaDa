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
#include <math/fuzzy.h>

namespace LaDa
{
  namespace CE
  {
    //! \brief Contains mapping policies for fixed-lattice separable functions.
    //! \todo Use boost::mpl.
    namespace Mapping
    {
      //! \brief Allows different types of mapping from confs to coef parameters.
      //! \details This mapping is equivalent to VectorPlus, eg (1,0..), (0,1,....),
      //!           and so on.
      template< size_t DIM > class VectorPlus 
      {
        public:
          //! A D dimensional mapping.
          const static size_t D = DIM;
          //! \brief Applies a functor with the correct coefficients for that
          //!        configuration element. 
          template< class T_OP, class T_CONF, class T_ITCOEF, class T_OUT >
            static void apply( const T_OP &_op, const T_CONF &_conf,
                               const T_ITCOEF &_coef, T_OUT &_out )
            { _op( _out, *( _coef + typename T_ITCOEF::difference_type( _conf ) ) ); }
          //! Applies function itself.
          template< class T_CONF, class T_ITCOEF, class T_OUT >
            static void apply( const T_CONF &_conf, const T_ITCOEF &_coef, 
                               T_OUT &_out )
            { 
              namespace bl = boost::lambda;
              apply( bl::_1 *= bl::_2, _conf, _coef, _out );
            }
          //! Returns the normalized, coef-less vector.
          template< class T_CONF, class T_OUT >
            static void add_tovec( const T_CONF &_conf, T_OUT &_out, 
                                   typename T_OUT::value_type _s  
                                       = typename T_OUT::value_type(1) )
            { _out[ size_t( _conf ) ] += typename T_OUT::value_type(_s); }
          //! Randomizes coefficients.
          template< class T_ITCOEF >
            static void randomize( T_ITCOEF &_coef,
                                   typename T_ITCOEF::value_type _range )
            {
              typedef typename T_ITCOEF :: value_type t_Type;
              *_coef = t_Type( opt::math::rng() - 0.5e0 ) * _range + t_Type(1);
              *(_coef + 1) = t_Type( opt::math::rng() - 0.5e0 ) * _range + t_Type(1);
              for( size_t i(2); i < D; ++i )
                *(_coef + i) = t_Type( opt::math::rng() - 0.5e0 ) * _range + t_Type(1);
            }

          //! Normalizes vector.
          template< class T_ITCOEF, class T_NORM >
            static void normalize( T_ITCOEF &_coef, T_NORM &_out )
            {
              namespace bl = boost::lambda;
              T_NORM norm(0);
              std::for_each( _coef, _coef + D, bl::var(norm) += bl::_1 * bl::_1 );
              if( math::is_zero( norm) )
              {
                std::fill( _coef, _coef + D, 1e0 );
                return;
              } 
              norm = std::sqrt( norm ) / std::sqrt( D );
              _out *= norm;
              norm = T_NORM(1) / norm;
              std::for_each( _coef, _coef + D, bl::_1 *= bl::constant(norm) );
            }
          static size_t norm( size_t ) { return 1; } 
      };
      //! \brief Allows different types of mapping from confs to coef parameters.
      //! \details This mapping is equivalent to VectorPlus, with one constant
      //!          vector, and all other vectors with a single non-zero component.
      template< size_t DIM > class VectorDiff
      {
        public:
          //! A D dimensional mapping.
          const static size_t D = DIM;
          //! \brief Applies a functor with the correct coefficients for that
          //!        configuration element. 
          template< class T_OP, class T_CONF, class T_ITCOEF, class T_OUT >
            static void apply( const T_OP &_op, const T_CONF &_conf,
                               const T_ITCOEF &_coef, 
                               T_OUT &_out )
            {
              typedef typename T_ITCOEF::difference_type t_size_t;
              typename T_ITCOEF::value_type coef( *_coef );
              if( not math::is_zero( _conf ) )
                coef += *( _coef + t_size_t( _conf ) );
              _op( _out, coef );
            }
          //! Applies functions with appropriate coef.
          template< class T_CONF, class T_ITCOEF, class T_OUT >
            static void apply( const T_CONF &_conf,
                               const T_ITCOEF &_coef, 
                               T_OUT &_out )
            {
              namespace bl = boost::lambda;
              apply( bl::_1 *= bl::_2, _conf, _coef, _out );
            }
          //! Returns the normalized, coef-less vector.
          template< class T_CONF, class T_OUT >
            static void add_tovec( const T_CONF &_conf, T_OUT &_out, 
                                   const typename T_OUT::value_type _s 
                                           = typename T_OUT::value_type(1) )
            {
              _out[0] += typename T_OUT::value_type(_s);
              if( math::is_zero( _conf ) ) return;
              _out[ size_t( _conf ) ] +=  typename T_OUT::value_type(_s); 
            }
          template< class T_ITCOEF >
            static void randomize( T_ITCOEF &_coef,
                                   typename T_ITCOEF::value_type _range )
            {
              typedef typename T_ITCOEF :: value_type t_Type;
              *_coef = t_Type( opt::math::rng() - 0.5e0 ) * _range + t_Type(1);
              *( _coef + 1 ) = t_Type( opt::math::rng() - 0.5e0 ) * _range;
              for( size_t i(2); i < D; ++i )
                *(_coef + i) = t_Type( opt::math::rng() - 0.5e0 ) * _range;
            }

          //! Normalizes vector.
          template< class T_ITCOEF, class T_NORM >
            static void normalize( T_ITCOEF &_coef,  T_NORM &_out )
            {
              namespace bl = boost::lambda;
              types::t_real norm( *_coef );
              if( math::is_zero( norm) ) return;
              norm = std::abs( norm );
              _out *= norm;
              norm = T_NORM(1) / norm;
              std::for_each( _coef, _coef + D, bl::_1 *= bl::constant(norm) );
            }
          static size_t norm( size_t _i ) { return _i == 0 ? 0: 1; } 
      };
    } // end of Mapping namespace.
  } // end of CE namespace.
} // namespace LaDa
#endif

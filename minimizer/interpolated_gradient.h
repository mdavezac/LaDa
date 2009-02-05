//
//  Version: $Id$
//
#ifndef _LADA_MINIMIZER_INTERPOLATED_GRADIENT_H_
#define _LADA_MINIMIZER_INTERPOLATED_GRADIENT_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <algorithm>
#include <cmath>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace LaDa
{
  namespace Minimizer
  {
    //! Interpolates gradient using set number of steps 2*\a _n + 1 and \a _stepsize.
    template< class T_FUNCTION, class T_CGS >
      void interpolated_gradient( const T_FUNCTION &_func,
                                  const typename T_FUNCTION :: t_Arg &_arg, 
                                  const T_CGS& _cgs,
                                  typename T_FUNCTION :: t_Return *const _grad,
                                  const size_t _n = 1,
                                  const typename T_FUNCTION :: t_Return _stepsize = 1e-3 )
      {
        if( _n == 0 ) 
        {
          gradient_fromdiff( _func, _arg, _grad, _stepsize );
          return;
        }
        namespace bnu = boost::numeric::ublas;
        typedef typename T_FUNCTION :: t_Arg t_Arg;
        typedef typename T_FUNCTION :: t_Return t_Return;
        typedef bnu::vector< t_Return > t_Vector;
        typedef bnu::matrix< t_Return > t_Matrix;
        const size_t nbgrad( _arg.size() );
        const size_t order(_n*2 + 1 );
        t_Arg arg( _arg );
        const t_Return gamma( _func( _arg ) );
        t_Vector b( order );
        t_Vector x( order,0 );
 
        // create A
        t_Matrix Aconst(order, order);
        for( size_t i(0); i < order; ++i )
          for( size_t j(0); j < order; ++j )
            Aconst(i,j) = std::pow( _stepsize * t_Return( int(i)-int(_n) ), j );

        const t_Matrix Atrans = bnu::trans( Aconst );
        Aconst = bnu::prec_prod( Atrans, Aconst  );
 
        for( size_t u(0); u < nbgrad; ++u )
        {
          const typename t_Arg :: value_type keep( _arg[u] );
          // create A and b
          for( size_t i(0); i < order; ++i ) 
            if( i == _n ) b(i) = gamma;
            else
            {
              arg[u] = keep + Atrans( 1, i );
              b(i) = _func( arg );
              arg[u] = keep;
            }
          t_Matrix A( Aconst );
          _cgs( A, x, bnu::prec_prod( Atrans, b ) );
          *(_grad + u ) = x[1];
        }
      }

    //! Computes gradient from two points only.
    template< class T_FUNCTION >
      void gradient_fromdiff( const T_FUNCTION &_func,
                              const typename T_FUNCTION :: t_Arg &_arg, 
                              typename T_FUNCTION :: t_Return *const _grad,
                              const typename T_FUNCTION :: t_Return _stepsize = 1e-3 )
      {
        namespace bnu = boost::numeric::ublas;
        typedef typename T_FUNCTION :: t_Arg t_Arg;
        typedef typename T_FUNCTION :: t_Return t_Return;
        typedef bnu::vector< t_Return > t_Vector;
        typedef bnu::matrix< t_Return > t_Matrix;
        const size_t nbgrad( _arg.size() );
        t_Arg arg( _arg );
        const t_Return gamma( _func( _arg ) );
        const t_Return invstep( t_Return(1) / _stepsize ); 

        for( size_t u(0); u < nbgrad; ++u )
        {
          const typename t_Arg :: value_type keep( _arg[u] );
          arg[u] += _stepsize;
          *(_grad + u ) = ( _func( arg ) - gamma ) * invstep;
          arg[u] = keep;
        }
      }

  } // namespace Minimizer

} // namespace LaDa
#endif

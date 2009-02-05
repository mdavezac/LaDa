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
                                  size_t _n = 1,
                                  const typename T_FUNCTION :: t_Return _stepsize = 1e-3 )
      {
        __DOASSERT( _n == 0, "Order too small.\n" )
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
        t_Vector x( order );
 
        // create A
        t_Matrix Aconst(order, order);
        for( size_t i(0); i < order; ++i )
          for( size_t j(0); j < order; ++j )
            if( j == 0 ) Aconst( i, j ) = 1;
            else if( i == _n ) Aconst( i,j ) = 0;
            else if( j % 2 and i < _n )
              Aconst(i,j) = -std::pow( _stepsize * t_Return( std::abs( int(i)-int(_n) ) ), j );
            else 
              Aconst(i,j) = std::pow( _stepsize * t_Return( std::abs( int(i)-int(_n) ) ), j );

        const t_Matrix Atrans = bnu::trans( Aconst );
        Aconst = bnu::prec_prod( Atrans, Aconst  );
 
        std::cout << "g: ";
        for( size_t i(0); i < nbgrad; ++i )
        {
          const typename t_Arg :: value_type keep( _arg[i] );
          std::cout << "keep: " << keep << "\n";
          std::fill( x.begin(), x.end(), 0 );
          // create A and b
          for( size_t j(0); j < order; ++j ) 
            if( j == _n ) b(j) == gamma;
            else
            {
              arg[i] = keep + Atrans( 1, j );
              b(j) = _func( arg );
              std::cout << "[ " << Atrans(1, j) << ", " <<  b(j) << " ]   ";
              arg[i] = keep;
            }
          t_Matrix A( Aconst );
          _cgs( A, x, bnu::prec_prod( Atrans, b ) );
          *(_grad + i ) = x[2];
          std::cout << x[2] << " ";
        }
        std::cout << "\n";
      }
  } // namespace Minimizer

} // namespace LaDa
#endif

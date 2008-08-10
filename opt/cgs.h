//
//  Version: $Id$
//
#ifndef _OPT_CGS_H_
#define _OPT_CGS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <utility>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "types.h"
#include "debug.h"
#include "fuzzy.h"

namespace Fitting
{
  //! \brief Conjugate Gradient Squared Method.
  //! \details Algorithm has been adapted from IML++ CGS to accept
  //!          boost::numberic::ublas array and matrices.
  class Cgs
  {
    //! Dummy preconditionner which does nothing.
    class DummyPrecond
    {
      public:
        //! The required interface for the preconditionner.
        template< class T_VECTOR >
          void solve( const T_VECTOR& _x ) const {}
    };

    public:
      //! Type of the return.
      typedef std::pair< types::t_real, types::t_unsigned > t_Return;
      //! Verbosity.
      bool verbose;
      //! Tolerance.
      types::t_real tolerance;
      //! Maximum number of iterations.
      types::t_unsigned itermax;

      //! Constructor.
      Cgs() : verbose( false ), tolerance( 1e-4 ), itermax(10) {}
      //! Copy Constructor.
      Cgs   ( const Cgs &_c )
          : verbose( _c.verbose ), tolerance( _c.tolerance ),
            itermax( _c.itermax ) {}
      //! Destructor.
      ~Cgs() {}

      //! \brief Performs a conjugate gradient squared minimization.
      //! \details Solves the unsymmetric linear system \a _A\a _x = \a _b
      //!          using the Conjugate Gradient Squared method.
      //!          Follows the algorithm described on p. 26 of the SIAM
      //!          Templates book.
      //! \tparam T_MATRIX should be a boost::numeric::matrix.
      //! \tparam T_VECTOR should be a boost::numeric::vector.
      //! \tparam T_PRECOND should be anything with a void solver( T_VECTOR& )
      //!                   member function. 
      //! \return a pair with the number of accomplished iterations in first,
      //!         and the convergence in second.
      template< class T_MATRIX, class T_VECTOR, class T_PRECOND >
      t_Return operator()( T_MATRIX &_A, T_VECTOR &_x, 
                           const T_VECTOR &_b, T_PRECOND &_precond ) const;
      
      //! \brief Performs a conjugate gradient squared minimization without
      //!        preconditionning.
      template< class T_MATRIX, class T_VECTOR >
      t_Return operator()( T_MATRIX &_A, T_VECTOR &_x, const T_VECTOR &_b ) const
        { DummyPrecond p; return operator()( _A, _x, _b, p ); }
    };
         

    template< class T_MATRIX, class T_VECTOR, class T_PRECOND >
      Cgs::t_Return Cgs :: operator()( T_MATRIX &_A, T_VECTOR &_x, 
                                       const T_VECTOR &_b, T_PRECOND &_precond ) const
      {
        __DEBUGTRYBEGIN
        __ASSERT( _A.size1() != _A.size2(), "matrix is not square.\n" )
        __ASSERT( _x.size() != _A.size2(), 
                     "Matrix and solution vector have different sizes.\n" 
                  << "  " << _x.size() << " != " << _A.size2() << "\n" )
        __ASSERT( _b.size() != _A.size1(), 
                     "Matrix and solution vector have different sizes.\n" 
                  << "  " << _b.size() << " != " << _A.size1() << "\n" )
        namespace ublas = boost::numeric::ublas;
        typedef T_VECTOR t_Vector;
        typedef typename t_Vector::value_type t_Type;
        t_Type resid, rho_1, rho_2, alpha, beta;
        t_Vector p, phat, q, qhat, vhat, u, uhat;

        if( verbose ) 
          std::cout << "Starting conjugate-gradient squared method: "
                    << "  Tolerance = " << tolerance
                    << " -- Max iter= " << itermax << "\n";
  
        types::t_real normb = ublas::norm_2( _b );
        t_Vector r = _b - ublas::prod( _A, _x);
        t_Vector rtilde = r;
  
        types::t_unsigned iter(0);
        if (normb == 0.0) normb = 1;
        
        resid = ublas::norm_2(r) / normb;
        if ( Fuzzy::leq( resid, tolerance) ) goto out;
  
        for (++iter; iter <= itermax; ++iter)
        {
          rho_1 = ublas::prec_inner_prod(rtilde, r);
          if (rho_1 == 0) 
          {
            resid = ublas::norm_2(r) / normb;
            break;
          }
          if (iter == 1) { u = r; p = u; }
          else
          {
            beta = rho_1 / rho_2;
            u = r + beta * q;
            p = u + beta * (q + beta * p);
          }
          phat = p; _precond.solve( phat );
          vhat = ublas::prod( _A, phat );
          alpha = rho_1 / ublas::prec_inner_prod(rtilde, vhat);
          q = u - alpha * vhat;
          uhat = u + q; _precond.solve( uhat );
          _x += alpha * uhat;
          qhat = ublas::prod( _A, uhat );
          r -= alpha * qhat;
          rho_2 = rho_1;
          resid = ublas::norm_2(r) / normb;
          if ( Fuzzy::leq( resid, tolerance) ) break;

          if( not verbose ) continue;

          std::cout << "  Iteration: " << iter 
                    << "  residual: " << resid << "\n";
        }
        
out:
        if( verbose ) 
          std::cout << "Cgs returns after " << iter
                    << (iter > 0 ? " iterations": " iteration")
                    << " and with a residual "  << resid << ".\n";
        return t_Return( resid, iter );
        __DEBUGTRYEND(, "Error encountered while solving Ax = b system.\n" )
      }

} // end of Fitting namespace.


#endif


//
//  Version: $Id$
//
#ifndef _GSL_LLSQ_H_
#define _GSL_LLSQ_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include <boost/lambda/lambda.hpp>

#include<opt/debug.h>
#include<opt/types.h>
#include<opt/gsl.h>
#include<tinyxml/tinyxml.h>

#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>

namespace LaDa
{
  //! Namespace for least-square fitting methods
  namespace Fitting
  {
    //! \cond
    namespace details
    {
      template< class T_FUNCTION >
        int gsl_f( const gsl_vector*, void*, gsl_vector* );
      template< class T_FUNCTION >
        int gsl_df( const gsl_vector*, void*, gsl_matrix* );
      template< class T_FUNCTION >
        int gsl_fdf( const gsl_vector* _x, void* _this, 
                     gsl_vector *_func, gsl_matrix* _jac );
    }
    //! \endcond
    // Non Linear least-square fitting using the gnu scientific libary.
    class NonLinearGsl 
    {
      public:
        //! Type of the solution vector.
        typedef std::vector< types::t_real > t_Vector;

      public:
        //! Required convergence.
        types::t_real tolerance;
        //! Verbosity
        bool verbose;
        //! Maximum number of iterations.
        types::t_unsigned itermax;


        //! Constructor.
        NonLinearGsl() : tolerance( 1e-4 ), verbose(false), itermax(0) {}
        //! Destructor.
        ~NonLinearGsl() {}

        //! Performs linear least-square fit.
        template< class T_FUNCTION >
        types::t_real operator()( T_FUNCTION &_func, t_Vector &_solution );

        //! Loads input from XML
        bool Load( const TiXmlElement &_node );
    };

    template< class T_FUNCTION >
      types::t_real NonLinearGsl::operator()( T_FUNCTION &_func, 
                                              t_Vector &_solution )
      {
        gsl_multifit_fdfsolver *solver = NULL;

        __DEBUGTRYBEGIN
        if( verbose )
          std::cout << " Starting Gsl Non-linear least-square fit.\n";

        int status;
        unsigned int i, iter = 0;
       
        gsl_multifit_function_fdf gsl_function;
        LaDa::Gsl::Vector x( _solution );
       
        gsl_function.f = &details::gsl_f<T_FUNCTION>;
        gsl_function.df = &details::gsl_df<T_FUNCTION>;
        gsl_function.fdf = &details::gsl_fdf<T_FUNCTION>;
        gsl_function.p = _solution.size();
        gsl_function.n = _func.dim();
        gsl_function.params = (void*) &_func;
       
       
        solver = gsl_multifit_fdfsolver_alloc( gsl_multifit_fdfsolver_lmsder, 
                                               _func.dim(), _solution.size() );
        gsl_multifit_fdfsolver_set (solver, &gsl_function, (gsl_vector*) x);
       
        do
        {
          ++iter;
          status = gsl_multifit_fdfsolver_iterate (solver);
       
          if (status)  break;
       
          status = gsl_multifit_test_delta ( solver->dx, solver->x, 
                                             tolerance, tolerance );
          if( verbose )
            std::cout << "  Iter " << iter << ": " 
                      << gsl_blas_dnrm2( solver->f ) << "\n";
        }
        while (    status == GSL_CONTINUE 
               and ( itermax == 0 or iter < itermax) );
       
       
        types::t_real result( gsl_blas_dnrm2( solver->f ) );
        if( verbose )
          std::cout << " Final Iteration " << iter << ": " << result << "\n";
        gsl_multifit_fdfsolver_free (solver);
        return result;

        __DEBUGTRYEND
        ( 
          if( solver) gsl_multifit_fdfsolver_free (solver);,
          "Error encountered while solving non-linear"
          " least-square fit with Gsl library.\n"
        )
      }


    //! \cond
    namespace details
    {
      template< class T_FUNCTION >
        int gsl_f( const gsl_vector* _x, void* _data, gsl_vector* _func )
        {
          T_FUNCTION *_this = (T_FUNCTION*) _data;
          (*_this)( _x->data, _func->data );
        }
      template< class T_FUNCTION >
        int gsl_df( const gsl_vector* _x, void* _data, gsl_matrix* _jac )
        {
          T_FUNCTION *_this = (T_FUNCTION*) _data;
          (*_this)( _x->data, _jac->data );
        }
      template< class T_FUNCTION >
        int gsl_fdf( const gsl_vector* _x, void* _data, 
                     gsl_vector *_func, gsl_matrix* _jac )
        { 
          gsl_f<T_FUNCTION>( _x, _data, _func );
          gsl_df<T_FUNCTION>( _x, _data, _jac );
          return GSL_SUCCESS;
        } 
    }
    //! \endcond

  }
} // namespace LaDa
#endif

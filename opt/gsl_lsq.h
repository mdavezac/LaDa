//
//  Version: $Id$
//
#ifndef _GSL_LLSQ_H_
#define _GSL_LLSQ_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>

#include<opt/debug.h>
#include <opt/types.h>
#include <tinyxml/tinyxml.h>

#include <gsl/gsl_block.h>
#include <gsl/gsl_multifit.h>

#ifdef _LADADEBUG
#include <boost/lambda/lambda.hpp>
#endif

namespace LaDa
{
  //! Namespace for least-square fitting methods
  namespace Fitting
  {
    // Linear least-square fitting using the gnu scientific libary.
    class Gsl 
    {
      public:
        //! Type of the Matrix/Linear Operator.
        typedef std::vector< types::t_real > t_Matrix;
        //! Type of the solution vector.
        typedef std::vector< types::t_real > t_Vector;

      public:
        //! Whether to include weights.
        bool doweights;
        //! Whether to perform single-value-decomposition flavor.
        bool dosvd;
        //! Required convergence.
        types::t_real tolerance;


        //! Constructor.
        Gsl() : doweights( false ), dosvd( false ), tolerance( 1e-4 ) {}
        //! Destructor.
        ~Gsl() {}
        //! Initialises matrix and vector.
        void init( t_Matrix &_mat, t_Vector &_v ) 
          { init_A(_mat); init_b( _v ); } 
        //! Initializes linear operator.
        void init_A( t_Matrix &_mat ) { A = _mat; }
        //! Initializes result vector.
        void init_b( t_Vector &_v ) { b = _v; }
        //! Initializes weight vector.
        void init_w( t_Vector &_w ) { w = _w; }

        //! Performs linear least-square fit.
        types::t_real operator()( t_Vector &_solution );

        //! Loads input from XML
        bool Load( const TiXmlElement &_node );

      public:
        //! Linear operator.
        t_Matrix A;
        //! Vector to which to fit.
        t_Vector b;
        //! Weight vector.
        t_Vector w; 
      
    };

  }
} // namespace LaDa
#endif

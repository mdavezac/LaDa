//
//  Version: $Id$
//
#ifndef _SEPARABLE_LLSQ_H_
#define _SEPARABLE_LLSQ_H_

#include <vector>

#include <opt/types.h>
#include <tinyxml/tinyxml.h>

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
      //! Whether to perform single-value-decomposition flavor.
      types::t_real tolerance;


      //! Constructor.
      Gsl() : doweight( false ), dosvd( false ), tolerance( 1e-4 ) {}
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

    protected:
      //! Linear operator.
      t_Matrix A;
      //! Vector to which to fit.
      t_Vector b;
      //! Weight vector.
      t_Vector w; 
    
  }

  namespace details
  {
    class GslVector
    {
      protected:
        gsl_bloc bloc;
        gsl_vector vector;

      public:
        GslVector( Gsl::t_Vector& _vec );
        ~GslVector() {};

        const gsl_vector* operator() { return gsl_vector; }
        const gsl_vector* const operator() const { return gsl_vector; }
    };
    class GslMatrix
    {
      protected:
        gsl_bloc bloc;
        gsl_Matrix matrix;

      public:
        GslVector( types::t_int _nrow, Gsl::t_Matrix& _vec );
        ~GslVector() {};

        const gsl_matrix* operator() { return &gsl_matrix; }
        const gsl_matrix* const operator() const { return &gsl_matrix; }
    };
  }
}
